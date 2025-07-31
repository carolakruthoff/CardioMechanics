

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
    CBSolverNewmarkBeta::Init(parameters,model);
    
    activeStress_ = new CBDataCtrl;
    activeStress_->Init(GetNumberOfElements());
    SetActiveStressDataSource(activeStress_);
    
    VecDuplicate(nodes_, &dx_);
    VecZeroEntries(dx_);
    contact_ = 0;
    nForces_= 0;
    
    // find Contact Handling Plugin
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
    
    
    mat_= parameters_->GetArray<TInt>("Solver.ActiveStressEstimator.Materials",{});
    std::set<TInt> nodesOfInterest;
    std::set<TInt> elementsOfInterest;
    
    for(auto& i : solidElements_)
    {
        if(dynamic_cast<CBElementSolidT10T4*>(i) != 0)
            dynamic_cast<CBElementSolidT10T4*>(i)->SwitchToT4();
        if(dynamic_cast<CBElementSolidT10RIT4*>(i) != 0)
            dynamic_cast<CBElementSolidT10RIT4*>(i)->SwitchToT4();
        
        if(mat_.size() != 0)
            if(find(mat_.begin(),mat_.end(),i->GetMaterialIndex())==mat_.end())
                continue;
        
        for(int j=0; j < i->GetNumberOfNodesIndices(); j++ )
        {
            nodesOfInterest.insert(i->GetNodeIndex(j));
        }
        elementsOfInterest.insert(i->GetIndex());
    }
    numNodesOfInterestIndices_ = 3*nodesOfInterest.size();
    numElementOfInterestIndices = elementsOfInterest.size();
    
    PetscInt* ni = new PetscInt[numNodesOfInterestIndices_];
    elementsOfInterestMapping_ = new PetscInt[numElementOfInterestIndices];
   
    
    std::set<TInt>::iterator it = nodesOfInterest.begin();
    
    PetscInt* nim = new PetscInt[numNodes_];
    
    for(int i=0; i < numNodes_; i++)
        nim[i] = -1;
    
    for(int i=0; i<nodesOfInterest.size(); i++)
    {
        nim[*it] = i;
        ni[3 * i + 0] = 3* *it + 0;
        ni[3 * i + 1] = 3* *it + 1;
        ni[3 * i + 2] = 3* *it + 2;
        it++;
    }
    
    TInt numTargets = parameters_->Get<TInt>("Solver.ActiveStressEstimator.NumberOfTargetNodes",-1);
    thikonov_ = parameters_->Get<TInt>("Solver.ActiveStressEstimator.Thikonov",0);
    if(thikonov_ != 0 && thikonov_ != 1 && thikonov_ != 2 )
        throw std::runtime_error("void CBSolverActiveStressEstimatorNewmarkBeta::Init(ParameterMap* parameters, CBModel* model");
    
    l1_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l1",1e-55);
    l2_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l2",1e-15);
    l3_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l3",1e-13);
    
    std::set<TInt> mn;
    
    srand (time(NULL));
    
    std::set<TInt> m = contact_->GetMasterNodesLocalIndices();
    std::set<TInt> m2;
    
    for(auto i=m.begin(); i != m.end(); i++)
        if(nodesOfInterest.find(*i)!=nodesOfInterest.end())
            m2.insert(*i);
    
    int r = m2.size() / numTargets;
    if(r <= 0 || numTargets == -1)
        mn = m2;
    else
    {
        for(auto i: m2)
        {
            if(rand() % r == 0)
                mn.insert(i);
        }
        
    }
    
    
    numMasterNodesIndices_ = 3*mn.size();
    masterNodesIndicesMapping_  = new PetscInt[numMasterNodesIndices_];
    masterNodesIndicesNodesOfInterestMapping_= new PetscInt[numMasterNodesIndices_];
    
    it = mn.begin();
    
    for(int i=0; i < mn.size(); i++)
    {
        PetscInt j = nim[*it];
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
    
    
    it = elementsOfInterest.begin();
    for(int i=0; i<elementsOfInterest.size(); i++)
    {
        elementsOfInterestMapping_[i] = *it;
        it++;
    }
    
    ISCreateGeneral(Petsc::Comm(), numMasterNodesIndices_, masterNodesIndicesMapping_, PETSC_COPY_VALUES, &masterNodesIndices_);
    ISCreateGeneral(Petsc::Comm(), numElementOfInterestIndices, elementsOfInterestMapping_, PETSC_COPY_VALUES, &elementsOfInterestIndices_);
    ISCreateGeneral(Petsc::Comm(), numNodesOfInterestIndices_, ni, PETSC_COPY_VALUES, &nodesOfInterestIndices_);
    VecDuplicate(nodes_, &dist_);
    for(auto& i : solidElements_)
    {
        if(dynamic_cast<CBElementSolidT10T4*>(i) != 0)
            dynamic_cast<CBElementSolidT10T4*>(i)->SwitchToT10();
        if(dynamic_cast<CBElementSolidT10RIT4*>(i) != 0)
            dynamic_cast<CBElementSolidT10RIT4*>(i)->SwitchToT10();
    }
    
    GenerateElementLaplacian();
    Petsc::CreateSeqVector(numElementOfInterestIndices, &lastSolution_);
    VecDuplicate(lastSolution_,&tau_);
    VecDuplicate(lastSolution_, &tmpSolution_);
    VecDuplicate(nodes_, &tmpNodes_);
    VecZeroEntries(lastSolution_);
    VecZeroEntries(tmpSolution_);
    VecZeroEntries(tau_);
    
    
    
//    bool* bc = GetNodesComponentsBoundaryConditionsGlobal();
//
//    for(int i=0; i < numNodesOfInterestIndices_; i++)
//        bc[ni[i]] = false;
//        
//    adapter_->LinkNodesComponentsBoundaryConditionsGlobal(bc);
    
    delete[] ni;
    
}

CBStatus CBSolverActiveStressEstimatorNewmarkBeta::EstimatorStep(PetscScalar time)
{
    for(auto& i : solidElements_)
    {
        if(dynamic_cast<CBElementSolidT10T4*>(i) != 0)
            dynamic_cast<CBElementSolidT10T4*>(i)->SwitchToT4();
        if(dynamic_cast<CBElementSolidT10RIT4*>(i) != 0)
            dynamic_cast<CBElementSolidT10RIT4*>(i)->SwitchToT4();
        
    }
    
    Vec subDist;
    
    Vec dtau;
    
    Mat dfdtau;
    Mat subdfdtau;
    
    UpdateGhostNodesAndLinkToAdapter();
    VecGetSubVector(dist_, masterNodesIndices_, &subDist);
    Mat subdfdx;
    
    Mat a;
    
    Petsc::CreateSeqVector(numElementOfInterestIndices, &dtau);
    Petsc::CreateSeqMatrix(3*numNodes_, GetNumberOfElements(), 5000, &dfdtau);
    
    adapter_->LinkNodalForcesActiveStressJacobian(dfdtau);
    formulation_->CalcNodalForcesActiveStressJacobian();
    MatAssemblyBegin(dfdtau, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(dfdtau, MAT_FINAL_ASSEMBLY);
    MatGetSubMatrix(dfdtau, nodesOfInterestIndices_, elementsOfInterestIndices_, MAT_INITIAL_MATRIX, &subdfdtau);
    
    if(staticEstimation_)
        CBSolver::CalcNodalForcesJacobian();
    else
    {
    
        TFloat ra, rb;
    
        ra = globalRayleighAlpha_;
        rb = globalRayleighBeta_;
    
        globalRayleighAlpha_ *= 0;
        globalRayleighBeta_ *= 0;
    
      auto timestep = timing_.GetTimeStep();
        VecAXPBYPCZ(tmpDisplacement_, timestep, timestep * timestep * 0.5 * (1 - 2 * beta_), 0, velocity_, acceleration_);
        VecAXPBYPCZ(tmpVelocity_, 1.0, (1 - gamma_) * timestep, 0, velocity_, acceleration_);
        VecCopy(tmpDisplacement_, displacement_);
        CalcDampingMatrix();
        CBSolverNewmarkBeta::CalcNodalForcesJacobian(displacement_, nodalForcesJacobian_);
    
        globalRayleighAlpha_ = ra;
        globalRayleighBeta_ = rb;
    }
    
    MatAssemblyBegin(nodalForcesJacobian_, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(nodalForcesJacobian_, MAT_FINAL_ASSEMBLY);

    MatAssemblyBegin(nodalForcesJacobian_, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(nodalForcesJacobian_, MAT_FINAL_ASSEMBLY);
    
    MatGetSubMatrix(nodalForcesJacobian_, nodesOfInterestIndices_, nodesOfInterestIndices_, MAT_INITIAL_MATRIX, &subdfdx);
    Mat subdfdxT;
    MatTranspose(subdfdx,MAT_INITIAL_MATRIX,&subdfdxT);
    
    IS iperm;
    IS perm;
    
    Mat f;
    MatGetOrdering(subdfdxT,MATORDERINGND,&perm,&iperm);
    MatGetFactor(subdfdxT,MATSOLVERMUMPS,MAT_FACTOR_LU,&f);
    MatLUFactorSymbolic(f,subdfdx,perm,iperm,0);
    MatLUFactorNumeric(f,subdfdx,0);
    
    Mat inv;
    MatCreateSeqDense(Petsc::Comm(), numMasterNodesIndices_, numNodesOfInterestIndices_,0,&inv);
    // MatDuplicate(b, MAT_DO_NOT_COPY_VALUES, &inv);
    std::cout << "Dim" << numMasterNodesIndices_ << " x " << numNodesOfInterestIndices_ << "\n";
    
    Vec ba;
    Vec e;
    Petsc::CreateSeqVector(numNodesOfInterestIndices_, &ba);
    VecDuplicate(ba, &e);
    PetscInt* ind = new PetscInt[numNodesOfInterestIndices_];
    
    for(int i=0; i < numNodesOfInterestIndices_ ; i++)
        ind[i]=i;
    
    for(int i=0; i < numMasterNodesIndices_ ; i++)
    {
        VecZeroEntries(ba);
        VecSetValue(ba,masterNodesIndicesNodesOfInterestMapping_[i],1,INSERT_VALUES);
        
        VecAssemblyBegin(ba);
        VecAssemblyEnd(ba);
        VecZeroEntries(e);
        MatSolve(f,ba, e);
        if(((i*100)/numMasterNodesIndices_+1)%10==0)
            Petsc::print << "\xd\t" << ((i*100)/numMasterNodesIndices_+1) <<"%                          " ;
        PetscScalar* vals;
        VecGetArray(e, &vals);
        
        MatSetValues(inv, 1, &i, numNodesOfInterestIndices_, ind, vals, INSERT_VALUES);
        VecRestoreArray(e, &vals);
    }
    MatAssemblyBegin(inv, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(inv, MAT_FINAL_ASSEMBLY);
    
    VecDestroy(&ba);
    VecDestroy(&e);
    delete ind;
    
    MatMatMult(inv, subdfdtau, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &a);
    
    Mat aTa;
    MatTransposeMatMult(a, a, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &aTa);
    Vec dd;
    VecDuplicate(dtau, &dd);
    
    if(thikonov_ == 0)
    {
        VecSet(dd, l1_);
        MatDiagonalSet(aTa, dd, ADD_ALL_VALUES);
    }
    else if(thikonov_ == 1 || thikonov_ == 2)
        MatAXPY(aTa, l3_, elementLaplacian_,DIFFERENT_NONZERO_PATTERN);
    else
        throw std::runtime_error("CBStatus CBSolverActiveStressEstimatorNewmarkBeta::EstimatorStep(PetscScalar time): Thikonov must be 0 [L2-Norm],1 [L2-Norm Laplace dtau] or 2 [L2-Norm Laplace tau]\n");
    
    Mat aT;
    Vec d;
    VecDuplicate(dtau, &d);
    
    MatTranspose(a,MAT_INITIAL_MATRIX,&aT);
    VecScale(subDist,1);
    MatMult(aT, subDist, d);
    VecAXPY(d,l2_, lastSolution_);
    Vec tt;
    VecDuplicate(d, &tt);
    
    if(thikonov_ == 2)
    {
        MatMult(elementLaplacian_, tau_,tt);
        VecAXPY(d, l3_, tt);
    }
    
    VecDestroy(&tt);
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
    KSPSolve(ksp,d, dtau);
    PetscScalar* t;
    VecGetArray(dtau,&t);
    VecCopy(lastSolution_, tmpSolution_);
    VecZeroEntries(lastSolution_);
    
    for(int i=0; i < numElementOfInterestIndices ; i++)
    {
        TFloat val = activeStress_->Get(lastTime_, elementsOfInterestMapping_[i]);
        if((val - t[i]) < 0)
        {
            activeStress_->Set(time,elementsOfInterestMapping_[i],0);
            VecSetValue(tau_, i, -val, ADD_ALL_VALUES);
            VecSetValue(lastSolution_,i,val,INSERT_ALL_VALUES);
            
        }
        else if((val - t[i]) > 2e5)
        {
            activeStress_->Set(time,elementsOfInterestMapping_[i],2e5);
            VecSetValue(tau_, i, -val + 2e5, ADD_ALL_VALUES);
            VecSetValue(lastSolution_,i,val+2e5,INSERT_ALL_VALUES);
        }
        else if(t[i] > 2000)
        {
            activeStress_->Set(time,elementsOfInterestMapping_[i],val - 2000);
            VecSetValue(tau_, i, -2000, ADD_ALL_VALUES);
            VecSetValue(lastSolution_,i,2000,INSERT_ALL_VALUES);
        }
        else if(t[i] < -2000)
        {
            activeStress_->Set(time,elementsOfInterestMapping_[i],val + 2000);
            VecSetValue(tau_, i, 2000, ADD_ALL_VALUES);
            VecSetValue(lastSolution_,i,-2000,INSERT_ALL_VALUES);
        }
        else
        {
            activeStress_->Set(time,elementsOfInterestMapping_[i],val - t[i]);
            VecSetValue(tau_, i, -t[i], ADD_ALL_VALUES);
            VecSetValue(lastSolution_,i,t[i],INSERT_ALL_VALUES);
        }
    }
    VecAssemblyBegin(lastSolution_);
    VecAssemblyEnd(lastSolution_);
    
    VecAssemblyBegin(tau_);
    VecAssemblyEnd(tau_);
    
    MatDestroy(&nodalForcesJacobian_);
    MatDestroy(&subdfdx);
    MatDestroy(&f);
    MatDestroy(&inv);
    MatDestroy(&a);
    MatDestroy(&aT);
    MatDestroy(&aTa);
    MatDestroy(&dfdtau);
    MatDestroy(&subdfdtau);
    MatDestroy(&subdfdxT);
    VecDestroy(&subDist);
    VecDestroy(&dtau);
    VecDestroy(&d);
    VecDestroy(&dd);
    ISDestroy(&iperm);
    ISDestroy(&perm);
    KSPDestroy(&ksp);
    
    for(auto& i : solidElements_)
    {
        if(dynamic_cast<CBElementSolidT10T4*>(i) != 0)
            dynamic_cast<CBElementSolidT10T4*>(i)->SwitchToT10();
        if(dynamic_cast<CBElementSolidT10RIT4*>(i) != 0)
            dynamic_cast<CBElementSolidT10RIT4*>(i)->SwitchToT10();
    }
    
    return CBStatus::SUCCESS;
}

CBStatus CBSolverActiveStressEstimatorNewmarkBeta::SolverStep(PetscScalar time, bool forceJacobianAndDampingRecalculation)
{
    Vec prevNodes;
    Vec target;
    Vec velo;
    Vec acc;
    Vec ds;
 
    for(auto& i : solidElements_)
    {
        if(dynamic_cast<CBElementSolidT10T4*>(i) != 0)
            dynamic_cast<CBElementSolidT10T4*>(i)->SwitchToT10();
        if(dynamic_cast<CBElementSolidT10RIT4*>(i) != 0)
            dynamic_cast<CBElementSolidT10RIT4*>(i)->SwitchToT10();
    }
    
    double alpha =  parameters_->Get<TFloat>("Solver.ActiveStressEstimator.Alpha",1e6);
    double beta = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.Beta",1e8);
    
    VecDuplicate(nodes_, &velo);
    VecDuplicate(nodes_, &acc);
    VecDuplicate(nodes_, &prevNodes);
    
    VecCopy(nodes_, prevNodes);
    VecCopy(velocity_, velo);
    VecCopy(acceleration_, acc);
    
    UpdateActiveStress(currentTime_);
    UpdateGhostNodesAndLinkToAdapter();
    CreateNodesJacobianAndLinkToAdapter();
    
    contact_->SwitchOn();
    contact_->SetAlpha(alpha);
    
    CBStatus rc;
    rc = CBSolverNewmarkBeta::SolverStep(time);
    if(rc != CBStatus::SUCCESS)
        return rc;
    
    if(time == timing_.GetStartTime())
        return rc;
    
    VecCopy(nodes_, tmpNodes_);
    
    UpdateGhostNodesAndLinkToAdapter();
    CreateNodesJacobianAndLinkToAdapter();
    
    contact_->SetAlpha(beta);
    contact_->SwitchOn();
    
    if(alpha != beta)
        rc = CBSolverNewmarkBeta::SolverStep(time);
    
    if(rc != CBStatus::SUCCESS)
        return rc;
    
    VecDuplicate(nodes_, &target);
    VecDuplicate(nodes_, &ds);
    VecCopy(nodes_,target);
    
    UpdateGhostNodesAndLinkToAdapter();
    CreateNodesJacobianAndLinkToAdapter();
    
    contact_->Apply(time); 
    contact_->GetMasterNodesDistancesToSlaveElements(&ds);

    VecCopy(prevNodes, nodes_);
    VecCopy(velo, velocity_);
    VecCopy(acc, acceleration_);
    
    UpdateGhostNodesAndLinkToAdapter();
        
    VecCopy(target, dist_);
    VecAXPY(dist_, -1, tmpNodes_);
    VecAXPY(dist_,1,ds);
    Vec subDist;
    PetscScalar norm=0;
    VecGetSubVector(dist_, masterNodesIndices_, &subDist);
    VecNorm(subDist, NORM_2, &norm);
    
    
    if(norm > 1e-15)
    {
        TFloat nu_ = 0.5;
        VecScale(tmpNodes_, (1-nu_));
        VecAXPY(tmpNodes_, nu_, prevNodes);
        VecCopy(tmpNodes_, nodes_);
        VecCopy(velo, velocity_);
        VecCopy(acc, acceleration_);
        contact_->SetAlpha(alpha);
        
        EstimatorStep(time);
        
        if(rc != CBStatus::SUCCESS)
            return rc;
        
    }
    contact_->SetAlpha(alpha);
    
    VecDestroy(&subDist);
    
    VecCopy(prevNodes, nodes_);
    VecCopy(velo, velocity_);
    VecCopy(acc, acceleration_);
    
    UpdateActiveStress(time);
    UpdateGhostNodesAndLinkToAdapter();
    CreateNodesJacobianAndLinkToAdapter();
    
    rc = CBSolverNewmarkBeta::SolverStep(time);
    
    VecDestroy(&target);
    VecDestroy(&prevNodes);
    
    if(rc == CBStatus::FAILED)
    {
        PetscScalar* t;
        VecGetArray(tmpSolution_,&t);
        for(int i=0; i < numElementOfInterestIndices ; i++)
            activeStress_->Set(time,elementsOfInterestMapping_[i],t[i]);
    }

    return rc;
    
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

void CBSolverActiveStressEstimatorNewmarkBeta::GenerateElementLaplacian()
{
    Mat cMat,subCMat;
    
    Vec diag;
    Vec b;
    
    Petsc::CreateSeqMatrix(GetNumberOfSolidElements(),GetNumberOfSolidElements(),100,&cMat);
    Petsc::CreateSeqVector(GetNumberOfSolidElements(), &diag);
    
    VecSet(diag,0);
    
    for(auto& e : solidElements_) {
        std::unordered_set<TInt> neighbors = adapter_->GetSolver()->GetModel()->GetSolidElementNeighborsCommonFace(e->GetIndex());
        for(auto& n : neighbors)
            MatSetValue(cMat,e->GetIndex(),(n),-1,INSERT_ALL_VALUES);
    }
    
    MatDiagonalSet(cMat, diag, INSERT_ALL_VALUES);
    VecDuplicate(diag, &b);
    VecSet(diag, 1);
    MatMult(cMat, diag, b);
    VecScale(b,-1);
    
    MatDiagonalSet(cMat, b, INSERT_ALL_VALUES);
    MatAssemblyBegin(cMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(cMat, MAT_FINAL_ASSEMBLY);
    MatGetSubMatrix(cMat, elementsOfInterestIndices_, elementsOfInterestIndices_, MAT_INITIAL_MATRIX, &subCMat);
    MatCreateSeqDense(Petsc::Comm(), numElementOfInterestIndices, numElementOfInterestIndices, 0,&elementLaplacian_);
    MatCopy(subCMat, elementLaplacian_, DIFFERENT_NONZERO_PATTERN);
    
    Mat tmp;
    MatTransposeMatMult(elementLaplacian_, elementLaplacian_,MAT_INITIAL_MATRIX , PETSC_DEFAULT, &tmp);
    MatCopy(tmp,elementLaplacian_,DIFFERENT_NONZERO_PATTERN);
    MatDestroy(&tmp);
    MatDestroy(&cMat);
    MatDestroy(&subCMat);
    
    VecDestroy(&diag);
    VecDestroy(&b);
    
}









