

//  CBSolverActiveStressEstimatorNewmarkBeta.cpp
//  CardioMechanics_unstable
//
//  Created by Thomas Fritz on 07.06.12.
//  Copyright (c) 2012 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
//

#include <stdlib.h>
#include <algorithm>
#include <iterator>

#include "CBSolverActiveStressEstimatorNewmarkBeta.h"
#include "CBContactHandling.h"
#include "CBElementSolidT10T4.h"
#include "CBElementSolidT10RIT4.h"

void CBSolverActiveStressEstimatorNewmarkBeta::Init(ParameterMap* parameters, CBModel* model)
{
    CBSolverNewmarkBeta::Init(parameters,model); // diff ActiveStressEstimator
    
    activeStress_ = new CBDataCtrl;
    activeStress_->Init(GetNumberOfElements());
    SetActiveStressDataSource(activeStress_);
    
    VecDuplicate(nodes_, &dx_);
    VecZeroEntries(dx_);
    // initialize objects of ContactHandling and DetermineNodalForces
    contact_ = 0;
    nForces_= 0;
    
    // find Contact Handling Plugin, throws error, if it's not there
    for(auto& p : plugins_)
    {
        if(dynamic_cast<CBContactHandling*>(p) != 0)
            contact_ = dynamic_cast<CBContactHandling*>(p);
        if(dynamic_cast<CBDetermineNodalForces*>(p) != 0)
            nForces_= dynamic_cast<CBDetermineNodalForces*>(p);
    }
    if(contact_ == 0)
        throw std::runtime_error("void CBParameterEstimator::Init(): ContactHandling plugin is needed for parameter estimation");
    //    if(nForces_ == 0)
    //        throw std::runtime_error("void CBParameterEstimator::Init(): DetermineNodalForces plugin is needed for parameter estimation");
    
    // mat_: list of material indices (consider only elements with these material)
    mat_= parameters_->GetArray<TInt>("Solver.ActiveStressEstimator.Materials",{});
    // allocate vectors
    std::set<TInt> nodesOfInterest;
    std::set<TInt> elementsOfInterest;

    // switch to element type T4, these are used for the Estimator Step
    for(auto& i : solidElements_)
    {
        if(dynamic_cast<CBElementSolidT10T4*>(i) != 0)
            dynamic_cast<CBElementSolidT10T4*>(i)->SwitchToT4();
        if(dynamic_cast<CBElementSolidT10RIT4*>(i) != 0)
            dynamic_cast<CBElementSolidT10RIT4*>(i)->SwitchToT4();
        // see if the elements have the materialnumber(s) specified in the list
        if(mat_.size() != 0) {
            if(find(mat_.begin(),mat_.end(),i->GetMaterialIndex())==mat_.end())
                continue;
        }
        // for each solid element of interest (of specified material) get all node indices of this element
        // save the node indices in nodesofInterest
        for(int j=0; j < i->GetNumberOfNodesIndices(); j++ )
        {
            nodesOfInterest.insert(i->GetNodeIndex(j));
        }
        elementsOfInterest.insert(i->GetIndex());
    }
    // there are 3 values for each node of interest -> n: amount of evaluation point components
    numNodesOfInterestIndices_ = 3*nodesOfInterest.size();
    numElementOfInterestIndices = elementsOfInterest.size();
    
    PetscInt* ni = new PetscInt[numNodesOfInterestIndices_];
    elementsOfInterestMapping_ = new PetscInt[numElementOfInterestIndices];
   
    
    std::set<TInt>::iterator it = nodesOfInterest.begin();
    // nim: number of nodes of interest mapping (?)
    PetscInt* nim = new PetscInt[numNodes_]; // diff CBSolverActiveStressEstimator PetscInt* nim instead of nim_
    
    for(int i=0; i < numNodes_; i++)
        nim[i] = -1;

    // For each node in nodesOfInterest, maps its global index to a local index in nim_
    // and fills ni with triplets of indices for x, y, z components
    for(int i=0; i<nodesOfInterest.size(); i++)
    {
        nim[*it] = i;
        ni[3 * i + 0] = 3* *it + 0;
        ni[3 * i + 1] = 3* *it + 1;
        ni[3 * i + 2] = 3* *it + 2;
        it++;
    }
    // retrieves number of target nodes if they are given in the .xml script
    TInt numTargets = parameters_->Get<TInt>("Solver.ActiveStressEstimator.NumberOfTargetNodes",-1);
    // diff to ActiveStressEstimator, get number (?) of thikonov parameters and throw error, if there are more than 2
    thikonov_ = parameters_->Get<TInt>("Solver.ActiveStressEstimator.Thikonov",0);
    if(thikonov_ != 0 && thikonov_ != 1 && thikonov_ != 2 )
        throw std::runtime_error("void CBSolverActiveStressEstimatorNewmarkBeta::Init(ParameterMap* parameters, CBModel* model");

    // li_ are penalty parameters for the thikonov regularization, they can all be set to 0
    // it is advised to set 1 to 2 of these parameters unequal to 0, but never all
    // these parameters unequal to zero are defined in the .xml script
    l1_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l1",1e-55);
    l2_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l2",1e-15);
    l3_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l3",1e-13);


    std::set<TInt> mn; // diff to ActiveStressEstimator

    // We want a deterministic initialization of the random number -> fixed seed for random number generation
    srand (time(NULL)); // diff to ActiveStressEstimator, time(NULL) instead of 0

    // m: amount of elements of interest
    std::set<TInt> m = contact_->GetMasterNodesLocalIndices();
    std::set<TInt> m2;

    // write only master nodes that are also nodes of interest in m2
    for(auto i=m.begin(); i != m.end(); i++)
        if(nodesOfInterest.find(*i)!=nodesOfInterest.end())
            m2.insert(*i);
    
    int r = m2.size() / numTargets;
    if(r <= 0 || numTargets == -1) // diff to ActiveStressEstimator there it is r<= 1
        mn = m2; // diff to ActiveStressEstimator there it is mn_ (<-global variable ?)
    else
    { // if numTargets is specified it randomly selects a subset of master nodes
        for(auto i: m2)
        {
            if(rand() % r == 0)
                mn.insert(i);
        }
        
    }
    
    // allocate arrays for mapping
    numMasterNodesIndices_ = 3*mn.size();
    masterNodesIndicesMapping_  = new PetscInt[numMasterNodesIndices_];
    masterNodesIndicesNodesOfInterestMapping_= new PetscInt[numMasterNodesIndices_];
    
    it = mn.begin();

    //Master nodes are mapped to their global and local indices
    for(int i=0; i < mn.size(); i++)
    {
        PetscInt j = nim[*it]; // diff ActiveStressEstimator, there it is nim_ (global PetscInt*)
        if(j == -1)
            throw std::runtime_error("Damn It");
        
        masterNodesIndicesMapping_[3 * i + 0] = 3* (*it) + 0;
        masterNodesIndicesMapping_[3 * i + 1] = 3* (*it) + 1;
        masterNodesIndicesMapping_[3 * i + 2] = 3* (*it) + 2;
        
        masterNodesIndicesNodesOfInterestMapping_[3 * i + 0] = 3* j + 0;
        masterNodesIndicesNodesOfInterestMapping_[3 * i + 1] = 3* j + 1;
        masterNodesIndicesNodesOfInterestMapping_[3 * i + 2] = 3* j + 2;
        it++;
    }
    
    // map elements of interest to their indices
    it = elementsOfInterest.begin();
    for(int i=0; i<elementsOfInterest.size(); i++)
    {
        elementsOfInterestMapping_[i] = *it;
        it++;
    }

    // Create PETSc index sets (IS) for master nodes, elements and nodes of interest (efficient data access)
    ISCreateGeneral(Petsc::Comm(), numMasterNodesIndices_, masterNodesIndicesMapping_, PETSC_COPY_VALUES, &masterNodesIndices_);
    ISCreateGeneral(Petsc::Comm(), numElementOfInterestIndices, elementsOfInterestMapping_, PETSC_COPY_VALUES, &elementsOfInterestIndices_);
    ISCreateGeneral(Petsc::Comm(), numNodesOfInterestIndices_, ni, PETSC_COPY_VALUES, &nodesOfInterestIndices_);
    // duplicates nodes_ into dist_ -> storing distances or residuals
    VecDuplicate(nodes_, &dist_);

    // Switch to T10 elements, this is usually done at the end of the estimator step
    // since the estimation is done with T4 elements, independent of the type of mesh elements
    for(auto& i : solidElements_)
    {
        if(dynamic_cast<CBElementSolidT10T4*>(i) != 0)
            dynamic_cast<CBElementSolidT10T4*>(i)->SwitchToT10();
        if(dynamic_cast<CBElementSolidT10RIT4*>(i) != 0)
            dynamic_cast<CBElementSolidT10RIT4*>(i)->SwitchToT10();
    }

    // set up laplacian matrix (spatial smoothing regularization)
    GenerateElementLaplacian();
    // creates and duplicates Petsc vectors
    Petsc::CreateSeqVector(numElementOfInterestIndices, &lastSolution_); // diff ActiveStressEstimator: lastSolution_ is ti_ there
    VecDuplicate(lastSolution_,&tau_);// diff ActiveStressEstimator: lastSolution_ is ti_ there, &tau_ is &ti1_ there
    VecDuplicate(lastSolution_, &tmpSolution_);// diff ActiveStressEstimator: lastSolution_ is ti_ there, &tmpSolution_ is &ti2_ there
    VecDuplicate(nodes_, &tmpNodes_);

    // diff ActiveStressEstimator next 3 lines do not appear there
    VecZeroEntries(lastSolution_);
    VecZeroEntries(tmpSolution_);
    VecZeroEntries(tau_);
    
    
    
//    bool* bc = GetNodesComponentsBoundaryConditionsGlobal();
//
//    for(int i=0; i < numNodesOfInterestIndices_; i++)
//        bc[ni[i]] = false;
//        
//    adapter_->LinkNodesComponentsBoundaryConditionsGlobal(bc);
    // deletes temporary array ni
    delete[] ni; // diff ActiveSTressEstimator, delete[] instead of delete
    
}

// diff ActiveStressEstimator: no UpdateMasterNodesOfInterest here

/// eki: the actual optimization / parameter estimation step, calls estimator + solver + processes the results
CBStatus CBSolverActiveStressEstimatorNewmarkBeta::SolverStep(PetscScalar time, bool forceJacobianAndDampingRecalculation)
{
    // vector allocation
    Vec prevNodes;
    Vec target;
    Vec velo; // diff ActiveStressEstimator
    Vec acc; // diff ActiveStressEstimator
    Vec ds;

    // switch element type to T10 (usually done after the estimator step)
    for (auto &i : solidElements_) {
        if (dynamic_cast<CBElementSolidT10T4 *>(i) != 0)
            dynamic_cast<CBElementSolidT10T4 *>(i)->SwitchToT10();
        if (dynamic_cast<CBElementSolidT10RIT4 *>(i) != 0)
            dynamic_cast<CBElementSolidT10RIT4 *>(i)->SwitchToT10();
    }

    // get alpha and beta from .xml-file
    double alpha =  parameters_->Get<TFloat>("Solver.ActiveStressEstimator.Alpha",1e6);
    double beta = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.Beta",1e8); // diff to ActiveStressEstimator

    // duplicate nodes_ into velo, acc and prevNodes
    VecDuplicate(nodes_, &velo); // diff ActiveStressEstimator
    VecDuplicate(nodes_, &acc); // diff ActiveStressEstimator
    VecDuplicate(nodes_, &prevNodes);

    VecCopy(nodes_, prevNodes);
    VecCopy(velocity_, velo); // diff ActiveStressEstimator
    VecCopy(acceleration_, acc); // diff ActiveStressEstimator

    // call UpdateActiveStress with the current timestep
    UpdateActiveStress(currentTime_);  // diff ActiveStressEstimator, currentTime_ instead of lastTime_
    UpdateGhostNodesAndLinkToAdapter();
    CreateNodesJacobianAndLinkToAdapter();

    // diff to ActiveSTressEstimator next 2 lines
    contact_->SwitchOn(); // switch on contact Problem
    contact_->SetAlpha(alpha); //set alpha

    CBStatus rc;
    // call parent SolverStep in CBSolverEquilibrium for the current timestep
    // the solution (of SolverStep) is the displacement to previous time step node coords
    // rc is the status of the whole process (SUCCESS vs FAILED)
    rc = CBSolverNewmarkBeta::SolverStep(time);
    if(rc != CBStatus::SUCCESS)
        return rc;

    if(time == timing_.GetStartTime())
        return rc;

    // diff ActiveStressEstimator, copy nodes_ into temporary nodes
    VecCopy(nodes_, tmpNodes_);
    // diff ActiveStressEstimator (next 4 lines)
    UpdateGhostNodesAndLinkToAdapter();
    CreateNodesJacobianAndLinkToAdapter();
    contact_->SetAlpha(beta); // set beta (or alpha to beta?)
    contact_->SwitchOn(); // switch on contact problem

    // diff ActiveStressEstimator (next 4 lines)
    if(alpha != beta)
        rc = CBSolverNewmarkBeta::SolverStep(time); //get status for SolverStep from NewmarkBeta for the current time, if alpha and beta are not equal
    if(rc != CBStatus::SUCCESS)
        return rc; //return status if failed

    // duplicate nodes into target and ds
    // create copies of current node positions for use in the estimation loop
    VecDuplicate(nodes_, &target);
    VecDuplicate(nodes_, &ds);
    VecCopy(nodes_,target);

    UpdateGhostNodesAndLinkToAdapter();
    CreateNodesJacobianAndLinkToAdapter();

    // diff ActiveStress Estimator, there is the following line here: VecZeroEntries(dist_);
    contact_->Apply(time);
    // diff ActiveStress Estimator, there is the following line here:  UpdateMasterNodesOfInterest();
    contact_->GetMasterNodesDistancesToSlaveElements(&ds); // diff ActiveStressEstimator &ds instead of &dist_

    // diff ActiveSTressEstimator (next 3 lines), copy these temporary variables into global (?) pointers
    VecCopy(prevNodes, nodes_);
    VecCopy(velo, velocity_);
    VecCopy(acc, acceleration_);

    UpdateGhostNodesAndLinkToAdapter(); // diff ActiveStressEstimator

    // diff ActiveStressEstimator next 3 lines
    VecCopy(target, dist_); //copys content of target into dist_
    VecAXPY(dist_, -1, tmpNodes_); // scaled vector addition dist_ = dist_ + (-1)*tmpNodes_
    VecAXPY(dist_,1,ds); // scaled vector addition dist_ = dist_ + 1*ds
    // calculate initial norm for stop criteria
    Vec subDist;
    PetscScalar norm=0;
    VecGetSubVector(dist_, masterNodesIndices_, &subDist);
    VecNorm(subDist, NORM_2, &norm);

    //diff ActiveStressEstimator: get absolute and relative tolerances here and start iterative estimation loop

    // diff ActiveStressEstimator: complete following loop
    if(norm > 1e-15)
    {
        TFloat nu_ = 0.5;
        VecScale(tmpNodes_, (1-nu_)); // scale tmpNodes: with (1-nu_)
        VecAXPY(tmpNodes_, nu_, prevNodes); // scaled vector addition tmpNodes_ = tmpNodes_ + nu_*prevNodes
        VecCopy(tmpNodes_, nodes_); // copy tmpNodes_ into nodes_
        VecCopy(velo, velocity_);
        VecCopy(acc, acceleration_);
        contact_->SetAlpha(alpha); // set alpha in ContactHandling

        EstimatorStep(time); // do estimator step

        if(rc != CBStatus::SUCCESS)
            return rc; // return status in case of a fail

    } // end estimation loop

    // diff ActiveStressEstimator: following 5 codelines
    contact_->SetAlpha(alpha); // diff ActiveStressEstimator: set alpha in ContactHandling

    VecDestroy(&subDist);

    VecCopy(prevNodes, nodes_);
    VecCopy(velo, velocity_);
    VecCopy(acc, acceleration_);


    UpdateActiveStress(time);
    UpdateGhostNodesAndLinkToAdapter();
    CreateNodesJacobianAndLinkToAdapter();

    rc = CBSolverNewmarkBeta::SolverStep(time); // diff ActiveStressEstimator, no solver step is done there

    // cleanup
    VecDestroy(&target);
    VecDestroy(&prevNodes);

    // Handle success and failure and update active stress data
    if(rc == CBStatus::FAILED)
    {
        PetscScalar* t;
        VecGetArray(tmpSolution_,&t);  // diff ActiveStressEstimator: tmpSolution_ instead of ti1_
        for(int i=0; i < numElementOfInterestIndices ; i++)
            activeStress_->Set(time,elementsOfInterestMapping_[i],t[i]);
        // diff ActiveStressEstimator: ti1_ is restored and copied, and in case of success the active stress is updated and ti1_ is copied and restored again
    }

    return rc;

}

CBStatus CBSolverActiveStressEstimatorNewmarkBeta::EstimatorStep(PetscScalar time) // diff ActiveStressEstimator: int step instead of EstimatorStep
{
    // switch to T4 elements (start of the active tension estimation)
    for(auto& i : solidElements_)
    {
        if(dynamic_cast<CBElementSolidT10T4*>(i) != 0)
            dynamic_cast<CBElementSolidT10T4*>(i)->SwitchToT4();
        if(dynamic_cast<CBElementSolidT10RIT4*>(i) != 0)
            dynamic_cast<CBElementSolidT10RIT4*>(i)->SwitchToT4();
        
    }

    // allocation of arrays
    Vec subDist;
    
    Vec dtau; // diff ActiveStressEstimator
    
    Mat dfdtau; // tangential stiffness matrix (change of nodal forces with respect to active stress)
    Mat subdfdtau; // submatrix of tangential stiffness matrix
    
    UpdateGhostNodesAndLinkToAdapter();
    // extract distances only for master nodes
    VecGetSubVector(dist_, masterNodesIndices_, &subDist);
    Mat subdfdx;  // diff ActiveStressEstimator
    
    Mat a; // product of reduced inverse of stiffness matrix and K_T describing change of all nodal forces
    
    Petsc::CreateSeqVector(numElementOfInterestIndices, &dtau);
    Petsc::CreateSeqMatrix(3*numNodes_, GetNumberOfElements(), 5000, &dfdtau);
    
    adapter_->LinkNodalForcesActiveStressJacobian(dfdtau);
    formulation_->CalcNodalForcesActiveStressJacobian();
    // formulation_ is an object of CBFormulation
    // stores tangential stiffness matrix in dfdtau
    MatAssemblyBegin(dfdtau, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(dfdtau, MAT_FINAL_ASSEMBLY);
    // extract submatrix of tangential stiffness matrix containing only nodes of interest
    MatGetSubMatrix(dfdtau, nodesOfInterestIndices_, elementsOfInterestIndices_, MAT_INITIAL_MATRIX, &subdfdtau);

    // diff ActiveStressEstimator: following if-else-statement
    if(staticEstimation_)
        CBSolver::CalcNodalForcesJacobian(); // call CalcNodalForcesJacobian, for static estimation
    else
    { // dynamic analysis
        // variable allocation
        TFloat ra, rb;
        // temporarily disable Raleigh damping (proportional damping) by saving current values ...
        ra = globalRayleighAlpha_;
        rb = globalRayleighBeta_;
        // ... and setting current values to zero
        globalRayleighAlpha_ *= 0;
        globalRayleighBeta_ *= 0;
    
        auto timestep = timing_.GetTimeStep(); // get current timestep
        VecAXPBYPCZ(tmpDisplacement_, timestep, timestep * timestep * 0.5 * (1 - 2 * beta_), 0, velocity_, acceleration_);
        // set tmpDisplacement = velocity_ * timestep + acceleration * timestep^2 * 0.5 *(1-2*beta)
        VecAXPBYPCZ(tmpVelocity_, 1.0, (1 - gamma_) * timestep, 0, velocity_, acceleration_);
        // set tmpVelocity_ = 1*velocity_ + (1 - gamma_)*timestep*acceleration + 0*tmpVelocity_
        VecCopy(tmpDisplacement_, displacement_); //update displacement_
        CalcDampingMatrix();
        CBSolverNewmarkBeta::CalcNodalForcesJacobian(displacement_, nodalForcesJacobian_); // CalcNodal forces with new displacement

        // reset damping variables
        globalRayleighAlpha_ = ra;
        globalRayleighBeta_ = rb;
    }

    // assemble nodalForcesJacobian_
    // Todo: does it have to be done twice?
    MatAssemblyBegin(nodalForcesJacobian_, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(nodalForcesJacobian_, MAT_FINAL_ASSEMBLY);

    MatAssemblyBegin(nodalForcesJacobian_, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(nodalForcesJacobian_, MAT_FINAL_ASSEMBLY);

    // get submatrix for nodesOfInterest: subdfdx
    MatGetSubMatrix(nodalForcesJacobian_, nodesOfInterestIndices_, nodesOfInterestIndices_, MAT_INITIAL_MATRIX, &subdfdx);
    Mat subdfdxT;
    // transpose subdfdx
    MatTranspose(subdfdx, MAT_INITIAL_MATRIX, &subdfdxT);

    // Compute LU factorization of subdfdxT using MUMPS
    IS iperm;
    IS perm;
    Mat f;
    MatGetOrdering(subdfdxT, MATORDERINGND, &perm, &iperm);
    MatGetFactor(subdfdxT, MATSOLVERMUMPS, MAT_FACTOR_LU, &f);
    MatLUFactorSymbolic(f, subdfdx, perm, iperm, 0); // prepare LU decomposition of subdfdx, result in f
    MatLUFactorNumeric(f, subdfdx, 0); // LU-Decomposition, invert

    // diff ActiveStressEstimator next 2 codelines -> appears before if-else statement in ActiveStressEstimator
    // dense matrix for inverse of reduced stiffness matrix
    Mat inv;
    MatCreateSeqDense(Petsc::Comm(), numMasterNodesIndices_, numNodesOfInterestIndices_,0,&inv);
    // MatDuplicate(b, MAT_DO_NOT_COPY_VALUES, &inv);
    std::cout << "Dim" << numMasterNodesIndices_ << " x " << numNodesOfInterestIndices_ << "\n";

    // Solves a series of linear equations to compute the inverse matrix
    Vec ba; // right hand side vector, set up for each master node, with 1 at corresponding position
    Vec e;
    Petsc::CreateSeqVector(numNodesOfInterestIndices_, &ba);
    VecDuplicate(ba, &e);
    PetscInt *ind = new PetscInt[numNodesOfInterestIndices_];
    
    for (int i = 0; i < numNodesOfInterestIndices_; i++)
        ind[i] = i;
    
    for (int i = 0; i < numMasterNodesIndices_; i++) {
        VecZeroEntries(ba);
        VecSetValue(ba, masterNodesIndicesNodesOfInterestMapping_[i],1,INSERT_VALUES); // get right hand side vector to calc inv
        
        VecAssemblyBegin(ba);
        VecAssemblyEnd(ba);
        VecZeroEntries(e);
        MatSolve(f, ba, e);
        // according to LeChat: solves subdfdxT * e = ba to get column of the inverse
        // inverse at this point is L
        if (((i*100)/numMasterNodesIndices_+1)%10 == 0)
            Petsc::print << "\xd\t" << ((i*100)/numMasterNodesIndices_+1) <<"%                          " ;
        PetscScalar *vals;
        VecGetArray(e, &vals);
        
        MatSetValues(inv, 1, &i, numNodesOfInterestIndices_, ind, vals, INSERT_VALUES);
        VecRestoreArray(e, &vals);
    }
    // assemble inverse matrix
    MatAssemblyBegin(inv, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(inv, MAT_FINAL_ASSEMBLY);

    // cleanup
    VecDestroy(&ba);
    VecDestroy(&e);
    // diff ActiveStressEstimator: even more pointers are destroyed here (&f, &subdfdxT, &iperm, &subdfdx, &perm)
    delete ind;

    // diff ActiveStressEstimator: tan operation reassembles the inv_, if inv_=0 and reuses the previous inverse if applicable

    // compute a = inv * subdfdtau = L * k^_d
    MatMatMult(inv, subdfdtau, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &a);
    
    Mat aTa;
    MatTransposeMatMult(a, a, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &aTa);


    // diff ActiveStressEstimator: d (right hand side vector) is initialized additionally with d = A^T g
    Vec dd;
    VecDuplicate(dtau, &dd);



    // diff ActiveStressEstimator: all following if-else statements)
    // Apply tikhonov regularization terms to stabilize solution -> check  for right configuration of regularization terms
    // each regularization term adds a weighted identity or Laplacian matrix to aTa and adjusts d
    if(thikonov_ == 0)
    {// add diagonal matrix dd scaled by l1_ to aTa -> 0. order thikonov (L2 norm penalty)
        VecSet(dd, l1_);
        MatDiagonalSet(aTa, dd, ADD_ALL_VALUES);  // aTa = aTa + dd
    }
    else if(thikonov_ == 1 || thikonov_ == 2)
        // adds scaled elementLapaplican_ to aTa -> 1. order thikonov (penalty to smoothness of solution)
        MatAXPY(aTa, l3_, elementLaplacian_,DIFFERENT_NONZERO_PATTERN);
        // aTa = l3_*elementLaplacian + aTa
    else
        throw std::runtime_error("CBStatus CBSolverActiveStressEstimatorNewmarkBeta::EstimatorStep(PetscScalar time): Thikonov must be 0 [L2-Norm],1 [L2-Norm Laplace dtau] or 2 [L2-Norm Laplace tau]\n");
    
    Mat aT;
    Vec d;
    VecDuplicate(dtau, &d);
    
    MatTranspose(a,MAT_INITIAL_MATRIX,&aT);
    VecScale(subDist,1);
    MatMult(aT, subDist, d);
    VecAXPY(d,l2_, lastSolution_);
    // d = aT * subDist  + l2_*lastSolution -> gradient of data misfit term plus regularization term
    Vec tt;
    VecDuplicate(d, &tt);
    
    if(thikonov_ == 2)
    { // add additional Laplacian term -> further penalizing smoothness of tau_
        MatMult(elementLaplacian_, tau_,tt); // tt = elementLaplacian_ * tau_
        VecAXPY(d, l3_, tt); // d = d + l3_*tt

    }
    // diff ActiveStressEstimator: cleanup
    VecDestroy(&tt);
    // Set up linear solver (KSP) with direct LU factorization
    // uses MUMPS
    KSP ksp;
    KSPCreate(Petsc::Comm(), &ksp);
    PC pc;
    KSPGetPC(ksp, &pc);
    KSPSetType(ksp, "preonly");
    PCFactorSetMatSolverPackage(pc,"mumps");
    KSPSetFromOptions(ksp);
    PCSetType(pc, PCLU);
    
    KSPSetOperators(ksp,aTa,aTa);
    KSPSetUp(ksp);
    // solves linear system aTa * dtau = d to find changes in active stress (dtau)
    KSPSolve(ksp,d, dtau);

    PetscScalar* t;
    // diff ActiveStressEstimator: PetscScalar *t2 is initialized here
    VecGetArray(dtau,&t);
    // diff ActiveStressEstimator next 2 code lines
    VecCopy(lastSolution_, tmpSolution_);
    VecZeroEntries(lastSolution_);

    // update active stress values for each element of interest for the lastTimeStep
    // diff ActiveStressEstimator: everything inside for loop
    for(int i=0; i < numElementOfInterestIndices ; i++)
    {
        TFloat val = activeStress_->Get(lastTime_, elementsOfInterestMapping_[i]); // get current activeStress value
        if((val - t[i]) < 0) // if difference between current active Stress and target is negative ...
        {
            activeStress_->Set(time,elementsOfInterestMapping_[i],0); // ... stress is clamped to 0
            VecSetValue(tau_, i, -val, ADD_ALL_VALUES);
            VecSetValue(lastSolution_,i,val,INSERT_ALL_VALUES);
            
        }
        else if((val - t[i]) > 2e5) // if difference exceeds 2e5 ...
        {
            activeStress_->Set(time,elementsOfInterestMapping_[i],2e5); // ... stress is clamped to 2e5
            VecSetValue(tau_, i, -val + 2e5, ADD_ALL_VALUES);
            VecSetValue(lastSolution_,i,val+2e5,INSERT_ALL_VALUES);
        }
        else if(t[i] > 2000) // if target stress exceeds 2000
        {
            activeStress_->Set(time,elementsOfInterestMapping_[i],val - 2000); // stress is reduced by 2000
            VecSetValue(tau_, i, -2000, ADD_ALL_VALUES);
            VecSetValue(lastSolution_,i,2000,INSERT_ALL_VALUES);
        }
        else if(t[i] < -2000) // if target stress is below -2000
        {
            activeStress_->Set(time,elementsOfInterestMapping_[i],val + 2000); // stress is increased by 2000
            VecSetValue(tau_, i, 2000, ADD_ALL_VALUES);
            VecSetValue(lastSolution_,i,-2000,INSERT_ALL_VALUES);
        }
        else //default case: stress is updated to val - t[i]
        {
            activeStress_->Set(time,elementsOfInterestMapping_[i],val - t[i]);
            VecSetValue(tau_, i, -t[i], ADD_ALL_VALUES);
            VecSetValue(lastSolution_,i,t[i],INSERT_ALL_VALUES);
        }
    }

    // diff ActiveStressEstimator: t2[i] = t[i], restore dtau, ti_ from &t and &t2

    // diff ActiveStressEstimator: next 4 code lines
    VecAssemblyBegin(lastSolution_);
    VecAssemblyEnd(lastSolution_); // combines all local contributions into global vector
    
    VecAssemblyBegin(tau_);
    VecAssemblyEnd(tau_); // combines all local contributions into global vector

    // cleanup
    MatDestroy(&nodalForcesJacobian_);
    MatDestroy(&subdfdx); // diff ActiveStressEstimator
    MatDestroy(&f); // diff ActiveStressEstimator
    MatDestroy(&inv);
    MatDestroy(&a);
    MatDestroy(&aT); // diff ActiveStressEstimator
    MatDestroy(&aTa);
    MatDestroy(&dfdtau);
    MatDestroy(&subdfdtau);
    MatDestroy(&subdfdxT); // diff ActiveStressEstimator
    VecDestroy(&subDist);
    VecDestroy(&dtau);
    VecDestroy(&d);
    VecDestroy(&dd);
    ISDestroy(&iperm); // diff ActiveStressEstimator
    ISDestroy(&perm); // diff ActiveStressEstimator
    KSPDestroy(&ksp);

    // switch to element type T10
    for(auto& i : solidElements_)
    {
        if(dynamic_cast<CBElementSolidT10T4*>(i) != 0)
            dynamic_cast<CBElementSolidT10T4*>(i)->SwitchToT10();
        if(dynamic_cast<CBElementSolidT10RIT4*>(i) != 0)
            dynamic_cast<CBElementSolidT10RIT4*>(i)->SwitchToT10();
    }
    
    return CBStatus::SUCCESS;
}



//CBStatus CBSolverActiveStressEstimatorNewmarkBeta::EstimatorStep(PetscScalar time)
//{
//    for(auto i=GetIteratorToSolidElementsBegin(); i != GetIteratorToSolidElementsEnd(); i++)
//    {
//        if(dynamic_cast<CBElementSolidT10T4*>(*i) != 0)
//            dynamic_cast<CBElementSolidT10T4*>(*i)->SwitchToT4();
//        if(dynamic_cast<CBElementSolidT10RIT4*>(*i) != 0)
//            dynamic_cast<CBElementSolidT10RIT4*>(*i)->SwitchToT4();
//        
//    }
//    
//    Vec subDist;
//    
//    Vec dtau;
//    
//    Mat dfdtau;
//    Mat subdfdtau;
//    
//    UpdateGhostNodesAndLinkToAdapter();
//    VecGetSubVector(dist_, masterNodesIndices_, &subDist);
//    Mat subdfdx;
//    
//    Mat a;
//    
//    Petsc::CreateSeqVector(numElementOfInterestIndices, &dtau);
//    Petsc::CreateSeqMatrix(3*numNodes_, GetNumberOfElements(), 5000, &dfdtau);
//    
//    adapter_->LinkNodalForcesActiveStressJacobian(dfdtau);
//    formulation_->CalcNodalForcesActiveStressJacobian();
//    MatAssemblyBegin(dfdtau, MAT_FINAL_ASSEMBLY);
//    MatAssemblyEnd(dfdtau, MAT_FINAL_ASSEMBLY);
//    MatGetSubMatrix(dfdtau, nodesOfInterestIndices_, elementsOfInterestIndices_, MAT_INITIAL_MATRIX, &subdfdtau);
//    Vec displacement;
//    
//    CBSolver::CalcNodalForcesJacobian();
//  
////    TFloat ra, rb;
////    
////    ra = globalRayleighAlpha_;
////    rb = globalRayleighBeta_;
////    
////    globalRayleighAlpha_ *= 0;
////    globalRayleighBeta_ *= 0;
////   
////    VecAXPBYPCZ(tmpDisplacement_, timeStep_, timeStep_ * timeStep_ * 0.5 * (1 - 2 * beta_), 0, velocity_, acceleration_);
////    VecAXPBYPCZ(tmpVelocity_, 1.0, (1 - gamma_) * timeStep_, 0, velocity_, acceleration_);
////    VecCopy(tmpDisplacement_, displacement_);
////    CalcDampingMatrix();
////    CBSolverNewmarkBeta::CalcNodalForcesJacobian(displacement_, nodalForcesJacobian_);
//// 
////    globalRayleighAlpha_ = ra;
////    globalRayleighBeta_ = rb;
//    
//    MatAssemblyBegin(nodalForcesJacobian_, MAT_FINAL_ASSEMBLY);
//    MatAssemblyEnd(nodalForcesJacobian_, MAT_FINAL_ASSEMBLY);
//    
//    MatGetSubMatrix(nodalForcesJacobian_, nodesOfInterestIndices_, nodesOfInterestIndices_, MAT_INITIAL_MATRIX, &subdfdx);
//    
//    Mat subdfdxT;
//    MatTranspose(subdfdx,MAT_INITIAL_MATRIX,&subdfdxT);
//    
//    IS iperm;
//    IS perm;
//    
//    Mat f;
//    MatGetOrdering(subdfdxT,MATORDERINGND,&perm,&iperm);
//    MatGetFactor(subdfdxT,MATSOLVERMUMPS,MAT_FACTOR_LU,&f);
//    MatLUFactorSymbolic(f,subdfdx,perm,iperm,0);
//    MatLUFactorNumeric(f,subdfdx,0);
//    
//    Mat inv;
//    MatCreateSeqDense(Petsc::Comm(), numMasterNodesIndices_, numNodesOfInterestIndices_,0,&inv);
//    // MatDuplicate(b, MAT_DO_NOT_COPY_VALUES, &inv);
//    std::cout << "Dim" << numMasterNodesIndices_ << " x " << numNodesOfInterestIndices_ << "\n";
//    
//    Vec ba;
//    Vec e;
//    Petsc::CreateSeqVector(numNodesOfInterestIndices_, &ba);
//    VecDuplicate(ba, &e);
//    PetscInt* ind = new PetscInt[numNodesOfInterestIndices_];
//    
//    for(int i=0; i < numNodesOfInterestIndices_ ; i++)
//        ind[i]=i;
//    
//    for(int i=0; i < numMasterNodesIndices_ ; i++)
//    {
//        VecZeroEntries(ba);
//        VecSetValue(ba,masterNodesIndicesNodesOfInterestMapping_[i],1,INSERT_VALUES);
//        
//        VecAssemblyBegin(ba);
//        VecAssemblyEnd(ba);
//        VecZeroEntries(e);
//        MatSolve(f,ba, e);
//        if(((i*100)/numMasterNodesIndices_+1)%10==0)
//            Petsc::print << "\xd\t" << ((i*100)/numMasterNodesIndices_+1) <<"%                          " ;
//        PetscScalar* vals;
//        VecGetArray(e, &vals);
//        
//        MatSetValues(inv, 1, &i, numNodesOfInterestIndices_, ind, vals, INSERT_VALUES);
//        VecRestoreArray(e, &vals);
//    }
//    MatAssemblyBegin(inv, MAT_FINAL_ASSEMBLY);
//    MatAssemblyEnd(inv, MAT_FINAL_ASSEMBLY);
//    
//    VecDestroy(&ba);
//    VecDestroy(&e);
//    delete ind;
//    
//    MatMatMult(inv, subdfdtau, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &a);
//    
//    Mat aTa;
//    MatTransposeMatMult(a, a, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &aTa);
//    Vec dd;
//    VecDuplicate(dtau, &dd);
//
//    if(thikonov_ == 0)
//    {
//        VecSet(dd, 1e-55);
//        MatDiagonalSet(aTa, dd, ADD_ALL_VALUES);
//    }
//    else if(thikonov_ == 1 || thikonov_ == 2)
//        MatAXPY(aTa, 1e-14, elementLaplacian_,DIFFERENT_NONZERO_PATTERN);
//    else
//        throw std::runtime_error("CBStatus CBSolverActiveStressEstimatorNewmarkBeta::EstimatorStep(PetscScalar time): Thikonov must be 0 [L2-Norm],1 [L2-Norm Laplace dtau] or 2 [L2-Norm Laplace tau]\n");
//    
//    Mat aT;
//    Vec d;
//    VecDuplicate(dtau, &d);
//    
//    MatTranspose(a,MAT_INITIAL_MATRIX,&aT);
//    VecScale(subDist,1);
//    MatMult(aT, subDist, d);
//    VecAXPY(d,1e-25 , lastSolution_);
//    Vec tt;
//    VecDuplicate(d, &tt);
//    
//    if(thikonov_ == 2)
//    {
//        MatMult(elementLaplacian_, tau_,tt);
//        VecAXPY(d, 1e-14, tt);
//    }
//    
//    VecDestroy(&tt);
//    KSP ksp;
//    KSPCreate(Petsc::Comm(), &ksp);
//    PC pc;
//    KSPGetPC(ksp, &pc);
//    KSPSetType(ksp, "preonly");
//    PCFactorSetMatSolverPackage(pc,"mumps");
//    KSPSetFromOptions(ksp);
//    PCSetType(pc, PCLU);
//    
//    KSPSetOperators(ksp,aTa,aTa,SAME_NONZERO_PATTERN);
//    KSPSetUp(ksp);
//    KSPSolve(ksp,d, dtau);
//    PetscScalar* t;
//    VecGetArray(dtau,&t);
//    VecCopy(lastSolution_, tmpSolution_);
//    VecZeroEntries(lastSolution_);
//    for(int i=0; i < numElementOfInterestIndices ; i++)
//    {
//        TFloat val = activeStress_->Get(lastTime_, elementsOfInterestMapping_[i]);
//        if((val - t[i]) < 0)
//        {
//            activeStress_->Set(time,elementsOfInterestMapping_[i],0);
//            VecSetValue(tau_, i, -val, ADD_ALL_VALUES);
//            VecSetValue(lastSolution_,i,val,INSERT_ALL_VALUES);
//            
//        }
//        else if((val - t[i]) > 2e5)
//        {
//            activeStress_->Set(time,elementsOfInterestMapping_[i],2e5);
//            VecSetValue(tau_, i, -val + 2e5, ADD_ALL_VALUES);
//            VecSetValue(lastSolution_,i,val+2e5,INSERT_ALL_VALUES);
//        }
//        else if(t[i] > 2000)
//        {
//            activeStress_->Set(time,elementsOfInterestMapping_[i],val - 2000);
//            VecSetValue(tau_, i, -2000, ADD_ALL_VALUES);
//            VecSetValue(lastSolution_,i,2000,INSERT_ALL_VALUES);
//        }
//        else if(t[i] < -2000)
//        {
//            activeStress_->Set(time,elementsOfInterestMapping_[i],val + 2000);
//            VecSetValue(tau_, i, 2000, ADD_ALL_VALUES);
//            VecSetValue(lastSolution_,i,-2000,INSERT_ALL_VALUES);
//        }
//        else
//        {
//            activeStress_->Set(time,elementsOfInterestMapping_[i],val - t[i]);
//            VecSetValue(tau_, i, -t[i], ADD_ALL_VALUES);
//            VecSetValue(lastSolution_,i,t[i],INSERT_ALL_VALUES);
//        }
//    }
//    VecAssemblyBegin(lastSolution_);
//    VecAssemblyEnd(lastSolution_);
//    
//    VecAssemblyBegin(tau_);
//    VecAssemblyEnd(tau_);
//    
//    MatDestroy(&nodalForcesJacobian_);
//    MatDestroy(&subdfdx);
//    MatDestroy(&f);
//    MatDestroy(&inv);
//    MatDestroy(&a);
//    MatDestroy(&aT);
//    MatDestroy(&aTa);
//    MatDestroy(&dfdtau);
//    MatDestroy(&subdfdtau);
//    MatDestroy(&subdfdxT);
//    VecDestroy(&subDist);
//    VecDestroy(&dtau);
//    VecDestroy(&d);
//    VecDestroy(&dd);
//    ISDestroy(&iperm);
//    ISDestroy(&perm);
//    KSPDestroy(&ksp);
//    for(auto i=GetIteratorToSolidElementsBegin(); i != GetIteratorToSolidElementsEnd(); i++)
//    {
//        if(dynamic_cast<CBElementSolidT10T4*>(*i) != 0)
//            dynamic_cast<CBElementSolidT10T4*>(*i)->SwitchToT10();
//        if(dynamic_cast<CBElementSolidT10RIT4*>(*i) != 0)
//            dynamic_cast<CBElementSolidT10RIT4*>(*i)->SwitchToT10();
//    }
//    
//}

/// eki: fills elementLaplacian_ with entries, actually just a connectivity matrix for regularization
void CBSolverActiveStressEstimatorNewmarkBeta::GenerateElementLaplacian()
{
    Mat cMat, subCMat; // C
    Vec diag;
    Vec b;
    // create sequential sparse matrix (cMat) with dimensions equal to number of solid elements
    Petsc::CreateSeqMatrix(GetNumberOfSolidElements(),GetNumberOfSolidElements(),100,&cMat);
    Petsc::CreateSeqVector(GetNumberOfSolidElements(), &diag);
    // create vector of same size as cMat, initialize to 0
    VecSet(diag,0);

    // iterate over all solid elements
    // build connectivity matrix where each element is connected to its neighbors
    for(auto& e : solidElements_) {
        // retrieve neighbor for each element (share common face)
        std::unordered_set<TInt> neighbors = adapter_->GetSolver()->GetModel()->GetSolidElementNeighborsCommonFace(e->GetIndex());
        // sets value of cMat at position (e->GetIndex(), n) to -1 for each neighbor n
        for(auto& n : neighbors)
            MatSetValue(cMat,e->GetIndex(),(n),-1,INSERT_ALL_VALUES);
    }
    // set diagonal entries of cMat to 0
    MatDiagonalSet(cMat, diag, INSERT_ALL_VALUES); // C = C + diag
    VecDuplicate(diag, &b); // b gets same size etc as diag, just allocates memory
    VecSet(diag, 1); // change all entries to 1
    // b = cMat*diag (row sums)
    MatMult(cMat, diag, b); // b = C * diag
    VecScale(b, -1); // b = -b
    // set diagonal entries of cMat to values in b
    MatDiagonalSet(cMat, b, INSERT_ALL_VALUES); // C = C + b
    // assemble cMat
    MatAssemblyBegin(cMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(cMat, MAT_FINAL_ASSEMBLY);
    // extract submatrix for elements of interest
    MatGetSubMatrix(cMat, elementsOfInterestIndices_, elementsOfInterestIndices_, MAT_INITIAL_MATRIX, &subCMat);
    MatCreateSeqDense(Petsc::Comm(), numElementOfInterestIndices, numElementOfInterestIndices, 0,&elementLaplacian_);
    // diff ActiveStressEstimator: MatDuplicate(subCMat, MAT_DO_NOT_COPY_VALUES, &elementLaplacian_);
    MatCopy(subCMat, elementLaplacian_, DIFFERENT_NONZERO_PATTERN);
    // compute Laplacian^T * Laplacian (often used in regularization or smoothing applications,
    // represents stiffness-like operator)
    Mat tmp; // diff ActiveStressEstimator
    MatTransposeMatMult(elementLaplacian_, elementLaplacian_, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &tmp);
    MatCopy(tmp,elementLaplacian_,DIFFERENT_NONZERO_PATTERN); // diff ActiveStressEstimator
    MatDestroy(&tmp); // diff ActiveStressEstimator
    // cleanup
    MatDestroy(&cMat);
    MatDestroy(&subCMat);
    VecDestroy(&diag);
    VecDestroy(&b);
}









