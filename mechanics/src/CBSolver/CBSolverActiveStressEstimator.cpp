/**@file CBSolverActiveStressEstimator.cpp
 * @brief <Brief (one-line) description here.>
 *
 * Please see the wiki for details on how to fill out this header:
 * https://intern.ibt.kit.edu/wiki/Document_a_IBT_C%2B%2B_tool
 *
 * @version 1.0.0
 *
 * @date Created <Your Name> (yyyy-mm-dd)
 *
 * @author Your Name\n
 *         Institute of Biomedical Engineering\n
 *         Karlsruhe Institute of Technology (KIT)\n
 *         http://www.ibt.kit.edu\n
 *         Copyright yyyy - All rights reserved.
 *
 * @see ...
 */

//  CBSolverActiveStressEstimator.cpp
//  CardioMechanics_unstable
//
//  Created by Thomas Fritz on 07.06.12.
//  Copyright (c) 2012 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
//

#include <stdlib.h>
#include <algorithm>
#include <iterator>

#ifdef  __APPLE__
# include <Accelerate/Accelerate.h>
#endif // ifdef  __APPLE__

#include "petscversion.h"
#ifndef MY_PETSC_VERSION
# define MY_PETSC_VERSION (PETSC_VERSION_MAJOR*10000 + PETSC_VERSION_MINOR*100+PETSC_VERSION_SUBMINOR)
#endif // ifndef MY_PETSC_VERSION
#if (MY_PETSC_VERSION >= 30600)
# include <petsc/private/matimpl.h> // for MatAXPY -- petsc >= 3.6
#else // if (MY_PETSC_VERSION >= 30600)
# include <petsc-private/matimpl.h> // for MatAXPY -- petsc <= 3.5
#endif // if (MY_PETSC_VERSION >= 30600)

#include "CBSolverActiveStressEstimator.h"
#include "CBContactHandling.h"
#include "CBElementSolidT10T4.h"
#include "CBElementSolidT10RIT4.h"

void CBSolverActiveStressEstimator::Init(ParameterMap *parameters, CBModel *model) {
  CBSolverEquilibrium::Init(parameters, model);

  // initialize active stress object
  activeStress_ = new CBDataCtrl;
  activeStress_->Init(GetNumberOfElements());
  SetActiveStressDataSource(activeStress_);

  //doubles nodes vector into dx and sets all values to zero -> probably dx used for displacement data
  VecDuplicate(nodes_, &dx_);
  VecZeroEntries(dx_);
  // initialize objects of ContactHandling and DetermineNodalForces
  contact_ = 0;
  nForces_ = 0;

  // find Contact Handling Plugin, throws error, if it's not there
  for (auto &p : plugins_) {
    if (dynamic_cast<CBContactHandling *>(p) != 0)
      contact_ = dynamic_cast<CBContactHandling *>(p);
    if (dynamic_cast<CBDetermineNodalForces *>(p) != 0)
      nForces_ = dynamic_cast<CBDetermineNodalForces *>(p);
  }
  if (contact_ == 0) {
    throw std::runtime_error(
            "void CBParameterEstimator::Init(): ContactHandling plugin is needed for parameter estimation");
  }

  //    if(nForces_ == 0)
  //        throw std::runtime_error("void CBParameterEstimator::Init(): DetermineNodalForces plugin is needed for parameter estimation");

  // mat_: list of material indices (consider only elements with these material)
  mat_ = parameters_->GetArray<TInt>("Solver.ActiveStressEstimator.Materials", {});
  // allocate vectors
  std::set<TInt> nodesOfInterest;
  std::set<TInt> elementsOfInterest;

  // switch to element type T4, these are used for the Estimator Step
  for (auto &i : solidElements_) {
    if (dynamic_cast<CBElementSolidT10T4 *>(i) != 0)
      dynamic_cast<CBElementSolidT10T4 *>(i)->SwitchToT4();
    if (dynamic_cast<CBElementSolidT10RIT4 *>(i) != 0)
      dynamic_cast<CBElementSolidT10RIT4 *>(i)->SwitchToT4();
  // see if the elements have the materialnumber(s) specified in the list
    if (mat_.size() != 0) {
      if (find(mat_.begin(), mat_.end(), i->GetMaterialIndex()) == mat_.end())
        continue;
    }
  // for each solid element of interest (of specified material) get all node indices of this element
  // save the node indices in nodesofInterest
    for (int j = 0; j < i->GetNumberOfNodesIndices(); j++) {
      nodesOfInterest.insert(i->GetNodeIndex(j));
    }
    elementsOfInterest.insert(i->GetIndex());
  }
  // there are 3 values for each node of interest -> n: amount of evaluation point components
  numNodesOfInterestIndices_ = 3*nodesOfInterest.size();
  numElementOfInterestIndices = elementsOfInterest.size();

  PetscInt *ni = new PetscInt[numNodesOfInterestIndices_];
  elementsOfInterestMapping_ = new PetscInt[numElementOfInterestIndices];


  std::set<TInt>::iterator it = nodesOfInterest.begin();

  // nim: number of nodes of interest mapping (?)
  nim_ = new PetscInt[numNodes_];

  for (int i = 0; i < numNodes_; i++)
    nim_[i] = -1;

  // For each node in nodesOfInterest, maps its global index to a local index in nim_
  // and fills ni with triplets of indices for x, y, z components
  for (int i = 0; i < nodesOfInterest.size(); i++) {
    nim_[*it] = i;
    ni[3 * i + 0] = 3* *it + 0;
    ni[3 * i + 1] = 3* *it + 1;
    ni[3 * i + 2] = 3* *it + 2;
    it++;
  }
  // retrieves number of target nodes if they are given in the .xml script
  TInt numTargets = parameters_->Get<TInt>("Solver.ActiveStressEstimator.NumberOfTargetNodes", -1);

  // li_ are penalty parameters for the thikonov regularization, they can all be set to 0
  // it is advised to set 1 to 2 of these parameters unequal to 0, but never all
  // these parameters unequal to zero are defined in the .xml script
  l1_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l1", 0);
  l2_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l2", 0);
  l3_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l3", 0);
  l4_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l4", 0);
  l5_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l5", 0);
  l6_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l6", 0);
  l7_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l7", 0);

  // We want a deterministic initialization of the random number -> fixed seed for random number generation
  srand(0);

  // m: amount of elements of interest
  std::set<TInt> m = contact_->GetMasterNodesLocalIndices();
  std::set<TInt> m2;

  // write only master nodes that are also nodes of interest in m2
  for (auto i = m.begin(); i != m.end(); i++)
    if (nodesOfInterest.find(*i) != nodesOfInterest.end())
      m2.insert(*i);

  int r = m2.size() / numTargets;
  if ((r <= 1) || (numTargets == -1)) {
    // mn: Masternodes of interest (?)
    mn_ = m2;
  } else { // if numTargets is specified it randomly selects a subset of master nodes
    for (auto i : m2) {
      if (rand() % r == 0)
        mn_.insert(i);
    }
  }

  // allocate arrays for mapping
  numMasterNodesIndices_ = 3*mn_.size();
  masterNodesIndicesMapping_  = new PetscInt[numMasterNodesIndices_];
  masterNodesIndicesNodesOfInterestMapping_ = new PetscInt[numMasterNodesIndices_];

  it = mn_.begin();

  //Master nodes are mapped to their global and local indices
  for (int i = 0; i < mn_.size(); i++) {
    PetscInt j = nim_[*it];
    if (j == -1)
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
  for (int i = 0; i < elementsOfInterest.size(); i++) {
    elementsOfInterestMapping_[i] = *it;
    it++;
  }

  // Create PETSc index sets (IS) for master nodes, elements and nodes of interest (efficient data access)
  ISCreateGeneral(
    Petsc::Comm(), numMasterNodesIndices_, masterNodesIndicesMapping_, PETSC_COPY_VALUES, &masterNodesIndices_);
  ISCreateGeneral(
    Petsc::Comm(), numElementOfInterestIndices, elementsOfInterestMapping_, PETSC_COPY_VALUES,
    &elementsOfInterestIndices_);
  ISCreateGeneral(Petsc::Comm(), numNodesOfInterestIndices_, ni, PETSC_COPY_VALUES, &nodesOfInterestIndices_);
  // duplicates nodes_ into dist_ -> storing distances or residuals
  VecDuplicate(nodes_, &dist_);

  // Switch to T10 elements, this is usually done at the end of the estimator step
  // since the estimation is done with T4 elements, independent of the type of mesh elements
  for (auto &i : solidElements_) {
    if (dynamic_cast<CBElementSolidT10T4 *>(i) != 0)
      dynamic_cast<CBElementSolidT10T4 *>(i)->SwitchToT10();
    if (dynamic_cast<CBElementSolidT10RIT4 *>(i) != 0)
      dynamic_cast<CBElementSolidT10RIT4 *>(i)->SwitchToT10();
  }

  // set up laplacian matrix (spatial smoothing regularization)
  GenerateElementLaplacian();
  // creates and duplicates Petsc vectors
  Petsc::CreateSeqVector(numElementOfInterestIndices, &ti_);
  VecDuplicate(ti_, &ti1_);
  VecDuplicate(ti_, &ti2_);
  VecDuplicate(nodes_, &tmpNodes_);

  //    bool* bc = GetNodesComponentsBoundaryConditionsGlobal();
  //
  //    for(int i=0; i < numNodesOfInterestIndices_; i++)
  //        bc[ni[i]] = false;
  //
  //    adapter_->LinkNodesComponentsBoundaryConditionsGlobal(bc);

  // deletes temporary array ni
  delete ni;
} // CBSolverActiveStressEstimator::Init

/// eki: take only nodes that possess a partner in the contact problem
void CBSolverActiveStressEstimator::UpdateMasterNodesOfInterest() {
  // update mapping for master nodes of interest

  // allocate arrays for master nodes
  std::set<TInt> mn;
  std::set<TInt> m = contact_->GetMasterNodesLocalIndices();
  std::set<TInt> m2;

  // write all master nodes, that are nodes of interest (in mn_) into the mn array
  for (auto i = m.begin(); i != m.end(); i++)
    if (mn_.find(*i) != mn_.end())
      mn.insert(*i);

  // array allocation for master node mapping
  numMasterNodesIndices_ = 3*mn.size();

  if (masterNodesIndices_)
    delete masterNodesIndicesMapping_;

  if (masterNodesIndicesNodesOfInterestMapping_)
    delete masterNodesIndicesNodesOfInterestMapping_;

  masterNodesIndicesMapping_  = new PetscInt[numMasterNodesIndices_];
  masterNodesIndicesNodesOfInterestMapping_ = new PetscInt[numMasterNodesIndices_];

  //Master nodes are mapped to their global and local indices
  auto it = mn.begin();

  for (int i = 0; i < mn.size(); i++) {
    PetscInt j = nim_[*it];
    if (j == -1)
      throw std::runtime_error("Damn It");

    masterNodesIndicesMapping_[3 * i + 0] = 3* (*it) + 0;
    masterNodesIndicesMapping_[3 * i + 1] = 3* (*it) + 1;
    masterNodesIndicesMapping_[3 * i + 2] = 3* (*it) + 2;

    masterNodesIndicesNodesOfInterestMapping_[3 * i + 0] = 3* j + 0;
    masterNodesIndicesNodesOfInterestMapping_[3 * i + 1] = 3* j + 1;
    masterNodesIndicesNodesOfInterestMapping_[3 * i + 2] = 3* j + 2;
    it++;
  }
} // CBSolverActiveStressEstimator::UpdateMasterNodesOfInterest

/// eki: the actual optimization / parameter estimation step, calls estimator + solver + processes the results
CBStatus CBSolverActiveStressEstimator::SolverStep(PetscScalar time, bool forceJacobianAndDampingRecalculation) {
  // vector allocation
  Vec prevNodes;
  Vec target;
  Vec ds;

  // switch element type to T10 (usually done after the estimator step)
  for (auto &i : solidElements_) {
    if (dynamic_cast<CBElementSolidT10T4 *>(i) != 0)
      dynamic_cast<CBElementSolidT10T4 *>(i)->SwitchToT10();
    if (dynamic_cast<CBElementSolidT10RIT4 *>(i) != 0)
      dynamic_cast<CBElementSolidT10RIT4 *>(i)->SwitchToT10();
  }

  // get alpha from .xml-file -> doesn't seem to be used
  double alpha =  parameters_->Get<TFloat>("Solver.ActiveStressEstimator.Alpha", 0);

  // duplicate nodes_ into prevNodes
  VecDuplicate(nodes_, &prevNodes);
  VecCopy(nodes_, prevNodes);

  // call UpdateActiveStress with the last timestep
  UpdateActiveStress(lastTime_);
  UpdateGhostNodesAndLinkToAdapter();
  CreateNodesJacobianAndLinkToAdapter();

  CBStatus rc;

  //    if(alpha != 0)
  //    {
  //        contact_->SwitchOn();
  //        contact_->SetAlpha(alpha);
  //
  //        rc = CBSolverEquilibrium::SolverStep(time);
  //        if(rc != CBStatus::SUCCESS)
  //            return rc;
  //
  //        if(time == startTime_)
  //            return rc;
  //    }
  //
  //    VecCopy(nodes_, tmpNodes_);
  //
  //    UpdateGhostNodesAndLinkToAdapter();
  //    CreateNodesJacobianAndLinkToAdapter();
  //
  //    contact_->SetAlpha(beta);
  //    contact_->SwitchOn();

  // call parent SolverStep in CBSolverEquilibrium for the current timestep
  // the solution (of SolverStep) is the displacement to previous time step node coords
  // rc is the status of the whole process (SUCCESS vs FAILED)
  rc = CBSolverEquilibrium::SolverStep(time);

  if (rc != CBStatus::SUCCESS)
    return rc;

  if (time == timing_.GetStartTime())
    return rc;

  //duplicate nodes into target and ds
  // create copies of current node positions for use in the estimation loop
  VecDuplicate(nodes_, &target);
  VecDuplicate(nodes_, &ds);
  VecCopy(nodes_, target);

  UpdateGhostNodesAndLinkToAdapter();
  CreateNodesJacobianAndLinkToAdapter();

  VecZeroEntries(dist_);

  contact_->Apply(time); // apply from ContactHandling
  UpdateMasterNodesOfInterest();
  contact_->GetMasterNodesDistancesToSlaveElements(&dist_);

  //    VecCopy(prevNodes, nodes_);
  //
  //    UpdateGhostNodesAndLinkToAdapter()t
  //
  //    VecCopy(target, dist_);
  //    VecAXPY(dist_, -1, tmpNodes_);
  //    VecAXPY(dist_,1,ds);
  // calculate initial norm for stop criteria
  Vec subDist;
  PetscScalar norm = 0;
  VecGetSubVector(dist_, masterNodesIndices_, &subDist);
  VecNorm(subDist, NORM_2, &norm);

  //    VecView(subDist, PETSC_VIEWER_STDOUT_SELF);
  Petsc::print << "Estimate ";
  Petsc::print << norm << "\n";
  int i = 0;
  // get absolute and relative tolerance from .xml-script for the stop criteria
  TFloat abs = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.AbsTol", 1e-4);
  TFloat rel = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.RelTol", 1e-4);

  // iterative estimation loop
  while (i < 5 && norm > abs) { // part of the stop criteria
    //        TFloat nu_ = 0.1;
    //        VecScale(tmpNodes_, (1-nu_));
    //        VecAXPY(tmpNodes_, nu_, prevNodes);
    //        VecCopy(tmpNodes_, nodes_);
    //        contact_->SetAlpha(alpha);

    // do estimator step
    EstimatorStep(time, i);

    //        contact_->SetAlpha(contact_->GetAlpha()/100);
    contact_->Apply(time); // apply from ContactHandling
    UpdateGhostNodesAndLinkToAdapter();
    CreateNodesJacobianAndLinkToAdapter();
    UpdateActiveStress(time);

    rc = CBSolverEquilibrium::SolverStep(time);

    UpdateGhostNodesAndLinkToAdapter();
    CreateNodesJacobianAndLinkToAdapter();

    VecZeroEntries(dist_);
    contact_->Apply(time); // apply from ContactHandling
    UpdateMasterNodesOfInterest();
    contact_->GetMasterNodesDistancesToSlaveElements(&dist_);
    VecGetSubVector(dist_, masterNodesIndices_, &subDist);
    TFloat lastNorm = norm;
    VecNorm(subDist, NORM_2, &norm);

    Petsc::print << norm << "\n";

    if (fabs(lastNorm - norm) < rel) // second part of stop criteria
      break;
    else
      i++;

    UpdateGhostNodesAndLinkToAdapter();
    CreateNodesJacobianAndLinkToAdapter();
  }

  //   VecDestroy(&subDist);
  //   VecCopy(prevNodes, nodes_);
  contact_->Apply(time);
  UpdateMasterNodesOfInterest();
  contact_->GetMasterNodesDistancesToSlaveElements(&dist_);

  //    VecCopy(prevNodes, nodes_);
  //
  //    UpdateGhostNodesAndLinkToAdapter()t
  //
  //    VecCopy(target, dist_);
  //    VecAXPY(dist_, -1, tmpNodes_);
  //    VecAXPY(dist_,1,ds);
  ExportNodesVectorData("Dist", dist_);

  UpdateActiveStress(time);
  UpdateGhostNodesAndLinkToAdapter();
  CreateNodesJacobianAndLinkToAdapter();

  //  rc = CBSolverEquilibrium::SolverStep(time);
  // cleanup
  VecDestroy(&target);
  VecDestroy(&prevNodes);

  // Handle success and failure and update active stress data
  if (rc == CBStatus::FAILED) {
    PetscScalar *t;
    VecGetArray(ti1_, &t);

    for (int i = 0; i < numElementOfInterestIndices; i++)
      activeStress_->Set(time, elementsOfInterestMapping_[i], -t[i]);

    VecRestoreArray(ti1_, &t);
    VecCopy(ti1_, ti_);
  } else {
    VecCopy(ti1_, ti2_);

    // fetch array from ti1_ as t
    PetscScalar *t;
    VecGetArray(ti1_, &t);

    for (int i = 0; i < numElementOfInterestIndices; i++)
      t[i] = -activeStress_->Get(time, elementsOfInterestMapping_[i]);
    VecRestoreArray(ti1_, &t);
    VecCopy(ti1_, ti_);
  }

  return rc;
} // CBSolverActiveStressEstimator::SolverStep

CBStatus CBSolverActiveStressEstimator::EstimatorStep(PetscScalar time, int step) {
  // switch to T4 elements (start of the active tension estimation)
  for (auto &i : solidElements_) {
    if (dynamic_cast<CBElementSolidT10T4 *>(i) != 0)
      dynamic_cast<CBElementSolidT10T4 *>(i)->SwitchToT4();
    if (dynamic_cast<CBElementSolidT10RIT4 *>(i) != 0)
      dynamic_cast<CBElementSolidT10RIT4 *>(i)->SwitchToT4();
  }

  // allocation of arrays
  Vec subDist;

  Mat dfdtau; // tangential stiffness matrix (change of nodal forces with respect to active stress)
  Mat subdfdtau; // submatrix of tangential stiffness matrix

  UpdateGhostNodesAndLinkToAdapter();
  // extract distances only for master nodes
  VecGetSubVector(dist_, masterNodesIndices_, &subDist);


  Mat a; // product of reduced inverse of stiffness matrix and K_T describing change of all nodal forces
  Vec dtau; // small change of active tension
  Petsc::CreateSeqVector(numElementOfInterestIndices, &dtau);
  Petsc::CreateSeqMatrix(3*numNodes_, GetNumberOfElements(), 5000, &dfdtau);

  adapter_->LinkNodalForcesActiveStressJacobian(dfdtau);
  formulation_->CalcNodalForcesActiveStressJacobian();
  // formulation_ is an object ob CBFormulation
  // stores tangential stiffness matrix in dfdtau
  MatAssemblyBegin(dfdtau, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(dfdtau, MAT_FINAL_ASSEMBLY);
  // extract submatrix of tangential stiffness matrix containing only nodes of interest
  MatGetSubMatrix(dfdtau, nodesOfInterestIndices_, elementsOfInterestIndices_, MAT_INITIAL_MATRIX, &subdfdtau);

  // dense matrix for inverse of reduced stiffness matrix
  Mat inv;
  MatCreateSeqDense(Petsc::Comm(), numMasterNodesIndices_, numNodesOfInterestIndices_, 0, &inv);

  // only enter this loop at step 0 or if the Inverse shouldn`t be reused
  if ((step == 0) || !(parameters_->Get<bool>("Solver.ActiveStressEstimator.ReUseInverse", false))) {
    Mat subdfdx; // subset of matrix describing change of all nodal forces with respect to nodal displacement
    CBSolver::CalcNodalForcesJacobian();
    MatAssemblyBegin(nodalForcesJacobian_, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(nodalForcesJacobian_, MAT_FINAL_ASSEMBLY);
    // calculation of subdfdx
    MatGetSubMatrix(nodalForcesJacobian_, nodesOfInterestIndices_, nodesOfInterestIndices_, MAT_INITIAL_MATRIX,
                    &subdfdx);
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

    // MatDuplicate(b, MAT_DO_NOT_COPY_VALUES, &inv);
    DCCtrl::print << "Dim " << numMasterNodesIndices_ << " x " << numNodesOfInterestIndices_ << "\n";

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
        DCCtrl::print << "\xd\t" << ((i*100)/numMasterNodesIndices_+1) <<"%                          ";
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
    MatDestroy(&f);
    MatDestroy(&subdfdxT);
    ISDestroy(&iperm);
    MatDestroy(&subdfdx);
    ISDestroy(&perm);
    delete[] ind;

    if (inv_ == 0)
      MatDuplicate(inv, MAT_DO_NOT_COPY_VALUES, &inv_);

    MatCopy(inv, inv_, DIFFERENT_NONZERO_PATTERN);
    MatAssemblyBegin(inv_, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(inv_, MAT_FINAL_ASSEMBLY);

  } else {// reuse previous inverse if applicable
    Petsc::debug << "Reusing previous matrix\n";
    MatCopy(inv_, inv, DIFFERENT_NONZERO_PATTERN);
    MatAssemblyBegin(inv, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(inv, MAT_FINAL_ASSEMBLY);
  }

  // compute a = inv * subdfdtau = L * k^_d
  MatMatMult(inv, subdfdtau, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &a);

  Mat aTa;
  MatTransposeMatMult(a, a, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &aTa);

  // Apply tikhonov regularization terms to stabilize solution
  // each regularization term adds a weighted identity or Laplacian matrix to aTa and adjusts d
  Vec d; // right-hand side vector
  VecDuplicate(dtau, &d);
  VecScale(subDist, 1);
  MatMultTranspose(a, subDist, d);  // ### d = A^T g, bzw. die rechte Seite um Ableitung 0 zu bestimmen

  Vec dd;
  Vec dd2; // diff ActiveStressTensorEstimator
  VecDuplicate(dtau, &dd);
  VecDuplicate(dtau, &dd2); // diff ActiveStressTensorEstimator

  // Applying regularization terms
  if (l1_ != 0) {
    VecSet(dd, l1_);
    MatDiagonalSet(aTa, dd, ADD_ALL_VALUES);  // aTa = aTa + dd
  }

  // eki: Was machen die Vektoren ti_, ti1, ti2 ? -> Daten aus letztem und vorletztem Zeitschritt
  if (l2_ != 0) {
    VecSet(dd, l2_);
    MatDiagonalSet(aTa, dd, ADD_ALL_VALUES);
    VecAXPY(d, -l2_, ti_);
  }

  auto timeStep = timing_.GetTimeStep();
  if (l3_ != 0) {
    VecSet(dd, l3_/ timeStep);
    MatDiagonalSet(aTa, dd, ADD_ALL_VALUES);
    VecAXPY(d, -l3_/timeStep, ti_);
    VecAXPY(d, l3_/timeStep, ti1_);
  }

  if (l4_ != 0) {
    VecSet(dd, l4_/ (timeStep*timeStep));
    MatDiagonalSet(aTa, dd, ADD_ALL_VALUES);

    VecAXPY(d, -l4_/(timeStep*timeStep), ti_);
    VecAXPY(d, 2*l4_/(timeStep*timeStep), ti1_);
    VecAXPY(d, -l4_/(timeStep*timeStep), ti2_);
  }


  if (l5_ != 0) {
    MatAXPY(aTa, l5_, lTl_, DIFFERENT_NONZERO_PATTERN);
  }

  if (l6_ != 0) {
    // sparse Matrix version
    MatAXPY(aTa, l6_, lTl_, DIFFERENT_NONZERO_PATTERN);  // aTa = aTa + l6_ * lTl_
    MatMult(lTl_, ti_, dd);                                    //  dd = lTl_ * ti_
    VecAXPY(d, -l6_, dd);                                     //   d = d - l6 * dd
  }

  if (l7_ != 0) {
    MatAXPY(aTa, l7_, lTl_, DIFFERENT_NONZERO_PATTERN);
    VecCopy(ti_, dd2);
    VecAXPY(dd2, -1, ti1_);
    MatMult(lTl_, dd2, dd);
    VecAXPY(d, -l7_, dd);
  }

  // -- -- -- -- -- -- -- -- -- --- -- --- ---
  // Set up linear solver (KSP) with direct LU factorization
  // uses MUMPS
  KSP ksp;
  KSPCreate(Petsc::Comm(), &ksp);
  PC pc;
  KSPGetPC(ksp, &pc);
  KSPSetType(ksp, "preonly");
  PCFactorSetMatSolverPackage(pc, "mumps");
  KSPSetFromOptions(ksp);
  PCSetType(pc, PCLU);

  KSPSetOperators(ksp, aTa, aTa);
  KSPSetUp(ksp);
  // solves linear system aTa * dtau = d to find changes in active stress (dtau)
  KSPSolve(ksp, d, dtau);

  // -- -- -- -- -- -- -- -- -- --- -- --- ---

  PetscScalar *t;
  PetscScalar *t2;

  VecGetArray(dtau, &t);
  VecGetArray(ti_, &t2);
  // update active stress values for each element of interest
  for (int i = 0; i < numElementOfInterestIndices; i++) {
    TFloat val = activeStress_->Get(lastTime_, elementsOfInterestMapping_[i]);

    if (t[i] > 10000)
      t[i] = 10000;
    else if (t[i] < -10000)
      t[i] = -10000;
    // applies bounds to change in stress to ensure numerical stability
    if ((val - t[i]) < 0) {
      t[i] = val;
      activeStress_->Set(time, elementsOfInterestMapping_[i], val - t[i]);
    } else if ((val - t[i]) > 2e5) {
      t[i] = val - 2e5;
      activeStress_->Set(time, elementsOfInterestMapping_[i], val - t[i]);
    }

    //        else if(dTau_[i] + t[i] > 5000)
    //        {
    //            t[i] = 5000 - dTau_[i];
    //            activeStress_->Set(time,elementsOfInterestMapping_[i], val - t[i]);
    //        }
    //
    //        else if(dTau_[i] + t[i] < -5000)
    //        {
    //            t[i] = -5000 + dTau_[i];
    //            activeStress_->Set(time,elementsOfInterestMapping_[i], val - t[i]);
    //        }
    else {
      activeStress_->Set(time, elementsOfInterestMapping_[i], val - t[i]);
    }
    t2[i] += t[i];

    //    std::cout << i << " " << t[i] << " " << t2[i] << "\n";
  }
  // update active stress data structure with new values
  VecRestoreArray(dtau, &t);
  VecRestoreArray(ti_, &t2);

  // cleanup
  MatDestroy(&nodalForcesJacobian_);
  MatDestroy(&inv);
  MatDestroy(&a);
  MatDestroy(&aTa);
  MatDestroy(&dfdtau);
  MatDestroy(&subdfdtau);
  VecDestroy(&subDist);
  VecDestroy(&dtau);
  VecDestroy(&d);
  VecDestroy(&dd);
  KSPDestroy(&ksp);
  // switch to element type T10
  for (auto &i : solidElements_) {
    if (dynamic_cast<CBElementSolidT10T4 *>(i) != 0)
      dynamic_cast<CBElementSolidT10T4 *>(i)->SwitchToT10();
    if (dynamic_cast<CBElementSolidT10RIT4 *>(i) != 0)
      dynamic_cast<CBElementSolidT10RIT4 *>(i)->SwitchToT10();
  }
  // update last estimator time step and return success status
  lastEstimatorTimeStep_ = time;
  lastTime_ = time;
  return CBStatus::SUCCESS;
} // CBSolverActiveStressEstimator::EstimatorStep

/// eki: fills elementLaplacian_ with entries, actually just a connectivity matrix for regularization
void CBSolverActiveStressEstimator::GenerateElementLaplacian() {
  Mat cMat, subCMat; // C
  Vec diag;
  Vec b;
  // create sequential sparse matrix (cMat) with dimensions equal to number of solid elements
  Petsc::CreateSeqMatrix(GetNumberOfSolidElements(), GetNumberOfSolidElements(), 100, &cMat);
  Petsc::CreateSeqVector(GetNumberOfSolidElements(), &diag);
  // create vector of same size as cMat, initialize to 0
  VecSet(diag, 0);

  // iterate over all solid elements
  // build connectivity matrix where each element is connected to its neighbors
  for (auto &e : solidElements_) {
    // retrieve neighbor for each element (share common face)
    std::unordered_set<TInt> neighbors = adapter_->GetSolver()->GetModel()->GetSolidElementNeighborsCommonFace(
      e->GetIndex());
    // sets value of cMat at position (e->GetIndex(), n) to -1 for each neighbor n
    for (auto &n : neighbors)
      MatSetValue(cMat, e->GetIndex(), n, -1, INSERT_ALL_VALUES);
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
  //extract submatrix for elements of interest
  MatGetSubMatrix(cMat, elementsOfInterestIndices_, elementsOfInterestIndices_, MAT_INITIAL_MATRIX, &subCMat);
  // create and duplicate element laplacian (copy subCMAt into it)
  MatCreateSeqAIJ(Petsc::Comm(), numElementOfInterestIndices, numElementOfInterestIndices, 0, 0, &elementLaplacian_);
  MatDuplicate(subCMat, MAT_DO_NOT_COPY_VALUES, &elementLaplacian_);
  MatCopy(subCMat, elementLaplacian_, DIFFERENT_NONZERO_PATTERN);
  // compute Laplacian^T * Laplacian (often used in regularization or smoothing applications,
  // represents stiffness-like operator)
  MatTransposeMatMult(elementLaplacian_, elementLaplacian_, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &lTl_);
  // cleanup
  MatDestroy(&cMat);
  MatDestroy(&subCMat);
  VecDestroy(&diag);
  VecDestroy(&b);
} // CBSolverActiveStressEstimator::GenerateElementLaplacian
