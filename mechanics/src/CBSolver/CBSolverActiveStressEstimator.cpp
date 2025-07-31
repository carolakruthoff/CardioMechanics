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

  activeStress_ = new CBDataCtrl;
  activeStress_->Init(GetNumberOfElements());
  SetActiveStressDataSource(activeStress_);

  VecDuplicate(nodes_, &dx_);
  VecZeroEntries(dx_);
  contact_ = 0;
  nForces_ = 0;

  // find Contact Handling Plugin
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


  mat_ = parameters_->GetArray<TInt>("Solver.ActiveStressEstimator.Materials", {});
  std::set<TInt> nodesOfInterest;
  std::set<TInt> elementsOfInterest;

  for (auto &i : solidElements_) {
    if (dynamic_cast<CBElementSolidT10T4 *>(i) != 0)
      dynamic_cast<CBElementSolidT10T4 *>(i)->SwitchToT4();
    if (dynamic_cast<CBElementSolidT10RIT4 *>(i) != 0)
      dynamic_cast<CBElementSolidT10RIT4 *>(i)->SwitchToT4();

    if (mat_.size() != 0) {
      if (find(mat_.begin(), mat_.end(), i->GetMaterialIndex()) == mat_.end())
        continue;
    }

    for (int j = 0; j < i->GetNumberOfNodesIndices(); j++) {
      nodesOfInterest.insert(i->GetNodeIndex(j));
    }
    elementsOfInterest.insert(i->GetIndex());
  }
  numNodesOfInterestIndices_ = 3*nodesOfInterest.size();
  numElementOfInterestIndices = elementsOfInterest.size();

  PetscInt *ni = new PetscInt[numNodesOfInterestIndices_];
  elementsOfInterestMapping_ = new PetscInt[numElementOfInterestIndices];


  std::set<TInt>::iterator it = nodesOfInterest.begin();

  nim_ = new PetscInt[numNodes_];

  for (int i = 0; i < numNodes_; i++)
    nim_[i] = -1;

  for (int i = 0; i < nodesOfInterest.size(); i++) {
    nim_[*it] = i;
    ni[3 * i + 0] = 3* *it + 0;
    ni[3 * i + 1] = 3* *it + 1;
    ni[3 * i + 2] = 3* *it + 2;
    it++;
  }

  TInt numTargets = parameters_->Get<TInt>("Solver.ActiveStressEstimator.NumberOfTargetNodes", -1);

  l1_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l1", 0);
  l2_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l2", 0);
  l3_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l3", 0);
  l4_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l4", 0);
  l5_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l5", 0);
  l6_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l6", 0);
  l7_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l7", 0);


  srand(0); // We want a deterministic initialization of the random number.

  std::set<TInt> m = contact_->GetMasterNodesLocalIndices();
  std::set<TInt> m2;

  for (auto i = m.begin(); i != m.end(); i++)
    if (nodesOfInterest.find(*i) != nodesOfInterest.end())
      m2.insert(*i);

  int r = m2.size() / numTargets;
  if ((r <= 1) || (numTargets == -1)) {
    mn_ = m2;
  } else {
    for (auto i : m2) {
      if (rand() % r == 0)
        mn_.insert(i);
    }
  }


  numMasterNodesIndices_ = 3*mn_.size();
  masterNodesIndicesMapping_  = new PetscInt[numMasterNodesIndices_];
  masterNodesIndicesNodesOfInterestMapping_ = new PetscInt[numMasterNodesIndices_];

  it = mn_.begin();

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


  it = elementsOfInterest.begin();
  for (int i = 0; i < elementsOfInterest.size(); i++) {
    elementsOfInterestMapping_[i] = *it;
    it++;
  }

  ISCreateGeneral(
    Petsc::Comm(), numMasterNodesIndices_, masterNodesIndicesMapping_, PETSC_COPY_VALUES, &masterNodesIndices_);
  ISCreateGeneral(
    Petsc::Comm(), numElementOfInterestIndices, elementsOfInterestMapping_, PETSC_COPY_VALUES,
    &elementsOfInterestIndices_);
  ISCreateGeneral(Petsc::Comm(), numNodesOfInterestIndices_, ni, PETSC_COPY_VALUES, &nodesOfInterestIndices_);
  VecDuplicate(nodes_, &dist_);
  for (auto &i : solidElements_) {
    if (dynamic_cast<CBElementSolidT10T4 *>(i) != 0)
      dynamic_cast<CBElementSolidT10T4 *>(i)->SwitchToT10();
    if (dynamic_cast<CBElementSolidT10RIT4 *>(i) != 0)
      dynamic_cast<CBElementSolidT10RIT4 *>(i)->SwitchToT10();
  }

  GenerateElementLaplacian();
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

  delete ni;
} // CBSolverActiveStressEstimator::Init

void CBSolverActiveStressEstimator::UpdateMasterNodesOfInterest() {
  std::set<TInt> mn;
  std::set<TInt> m = contact_->GetMasterNodesLocalIndices();
  std::set<TInt> m2;

  for (auto i = m.begin(); i != m.end(); i++)
    if (mn_.find(*i) != mn_.end())
      mn.insert(*i);


  numMasterNodesIndices_ = 3*mn.size();

  if (masterNodesIndices_)
    delete masterNodesIndicesMapping_;

  if (masterNodesIndicesNodesOfInterestMapping_)
    delete masterNodesIndicesNodesOfInterestMapping_;

  masterNodesIndicesMapping_  = new PetscInt[numMasterNodesIndices_];
  masterNodesIndicesNodesOfInterestMapping_ = new PetscInt[numMasterNodesIndices_];

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

CBStatus CBSolverActiveStressEstimator::SolverStep(PetscScalar time) {
  Vec prevNodes;
  Vec target;
  Vec ds;

  for (auto &i : solidElements_) {
    if (dynamic_cast<CBElementSolidT10T4 *>(i) != 0)
      dynamic_cast<CBElementSolidT10T4 *>(i)->SwitchToT10();
    if (dynamic_cast<CBElementSolidT10RIT4 *>(i) != 0)
      dynamic_cast<CBElementSolidT10RIT4 *>(i)->SwitchToT10();
  }

  double alpha =  parameters_->Get<TFloat>("Solver.ActiveStressEstimator.Alpha", 0);


  VecDuplicate(nodes_, &prevNodes);
  VecCopy(nodes_, prevNodes);

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
  //
  rc = CBSolverEquilibrium::SolverStep(time);

  if (rc != CBStatus::SUCCESS)
    return rc;

  if (time == timing_.GetStartTime())
    return rc;

  VecDuplicate(nodes_, &target);
  VecDuplicate(nodes_, &ds);
  VecCopy(nodes_, target);

  UpdateGhostNodesAndLinkToAdapter();
  CreateNodesJacobianAndLinkToAdapter();

  VecZeroEntries(dist_);

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
  Vec subDist;
  PetscScalar norm = 0;
  VecGetSubVector(dist_, masterNodesIndices_, &subDist);
  VecNorm(subDist, NORM_2, &norm);

  //    VecView(subDist, PETSC_VIEWER_STDOUT_SELF);
  Petsc::print << "Estimate ";
  Petsc::print << norm << "\n";
  int i = 0;
  TFloat abs = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.AbsTol", 1e-4);
  TFloat rel = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.RelTol", 1e-4);

  while (i < 5 && norm > abs) {
    //        TFloat nu_ = 0.1;
    //        VecScale(tmpNodes_, (1-nu_));
    //        VecAXPY(tmpNodes_, nu_, prevNodes);
    //        VecCopy(tmpNodes_, nodes_);
    //        contact_->SetAlpha(alpha);
    EstimatorStep(time, i);

    //        contact_->SetAlpha(contact_->GetAlpha()/100);
    contact_->Apply(time);
    UpdateGhostNodesAndLinkToAdapter();
    CreateNodesJacobianAndLinkToAdapter();
    UpdateActiveStress(time);

    rc = CBSolverEquilibrium::SolverStep(time);

    UpdateGhostNodesAndLinkToAdapter();
    CreateNodesJacobianAndLinkToAdapter();

    VecZeroEntries(dist_);
    contact_->Apply(time);
    UpdateMasterNodesOfInterest();
    contact_->GetMasterNodesDistancesToSlaveElements(&dist_);
    VecGetSubVector(dist_, masterNodesIndices_, &subDist);
    TFloat lastNorm = norm;
    VecNorm(subDist, NORM_2, &norm);

    Petsc::print << norm << "\n";

    if (fabs(lastNorm - norm) < rel)
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

  VecDestroy(&target);
  VecDestroy(&prevNodes);

  if (rc == CBStatus::FAILED) {
    PetscScalar *t;
    VecGetArray(ti1_, &t);

    for (int i = 0; i < numElementOfInterestIndices; i++)
      activeStress_->Set(time, elementsOfInterestMapping_[i], -t[i]);

    VecRestoreArray(ti1_, &t);
    VecCopy(ti1_, ti_);
  } else {
    VecCopy(ti1_, ti2_);
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
  for (auto &i : solidElements_) {
    if (dynamic_cast<CBElementSolidT10T4 *>(i) != 0)
      dynamic_cast<CBElementSolidT10T4 *>(i)->SwitchToT4();
    if (dynamic_cast<CBElementSolidT10RIT4 *>(i) != 0)
      dynamic_cast<CBElementSolidT10RIT4 *>(i)->SwitchToT4();
  }

  Vec subDist;


  Mat dfdtau;
  Mat subdfdtau;

  UpdateGhostNodesAndLinkToAdapter();
  VecGetSubVector(dist_, masterNodesIndices_, &subDist);


  Mat a;
  Vec dtau;
  Petsc::CreateSeqVector(numElementOfInterestIndices, &dtau);
  Petsc::CreateSeqMatrix(3*numNodes_, GetNumberOfElements(), 5000, &dfdtau);

  adapter_->LinkNodalForcesActiveStressJacobian(dfdtau);
  formulation_->CalcNodalForcesActiveStressJacobian();
  MatAssemblyBegin(dfdtau, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(dfdtau, MAT_FINAL_ASSEMBLY);
  MatGetSubMatrix(dfdtau, nodesOfInterestIndices_, elementsOfInterestIndices_, MAT_INITIAL_MATRIX, &subdfdtau);

  Mat inv;
  MatCreateSeqDense(Petsc::Comm(), numMasterNodesIndices_, numNodesOfInterestIndices_, 0, &inv);

  if ((step == 0) || !(parameters_->Get<bool>("Solver.ActiveStressEstimator.ReUseInverse", false))) {
    Mat subdfdx;
    CBSolver::CalcNodalForcesJacobian();
    MatAssemblyBegin(nodalForcesJacobian_, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(nodalForcesJacobian_, MAT_FINAL_ASSEMBLY);

    MatGetSubMatrix(nodalForcesJacobian_, nodesOfInterestIndices_, nodesOfInterestIndices_, MAT_INITIAL_MATRIX,
                    &subdfdx);
    Mat subdfdxT;
    MatTranspose(subdfdx, MAT_INITIAL_MATRIX, &subdfdxT);

    IS iperm;
    IS perm;

    Mat f;
    MatGetOrdering(subdfdxT, MATORDERINGND, &perm, &iperm);
    MatGetFactor(subdfdxT, MATSOLVERMUMPS, MAT_FACTOR_LU, &f);
    MatLUFactorSymbolic(f, subdfdx, perm, iperm, 0);
    MatLUFactorNumeric(f, subdfdx, 0);

    // MatDuplicate(b, MAT_DO_NOT_COPY_VALUES, &inv);
    DCCtrl::print << "Dim " << numMasterNodesIndices_ << " x " << numNodesOfInterestIndices_ << "\n";

    Vec ba;
    Vec e;
    Petsc::CreateSeqVector(numNodesOfInterestIndices_, &ba);
    VecDuplicate(ba, &e);
    PetscInt *ind = new PetscInt[numNodesOfInterestIndices_];

    for (int i = 0; i < numNodesOfInterestIndices_; i++)
      ind[i] = i;

    for (int i = 0; i < numMasterNodesIndices_; i++) {
      VecZeroEntries(ba);
      VecSetValue(ba, masterNodesIndicesNodesOfInterestMapping_[i], 1, INSERT_VALUES);

      VecAssemblyBegin(ba);
      VecAssemblyEnd(ba);
      VecZeroEntries(e);
      MatSolve(f, ba, e);
      if (((i*100)/numMasterNodesIndices_+1)%10 == 0)
        DCCtrl::print << "\xd\t" << ((i*100)/numMasterNodesIndices_+1) <<"%                          ";
      PetscScalar *vals;
      VecGetArray(e, &vals);

      MatSetValues(inv, 1, &i, numNodesOfInterestIndices_, ind, vals, INSERT_VALUES);
      VecRestoreArray(e, &vals);
    }
    MatAssemblyBegin(inv, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(inv, MAT_FINAL_ASSEMBLY);

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
  } else {
    Petsc::debug << "Reusing previous matrix\n";
    MatCopy(inv_, inv, DIFFERENT_NONZERO_PATTERN);
    MatAssemblyBegin(inv, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(inv, MAT_FINAL_ASSEMBLY);
  }


  MatMatMult(inv, subdfdtau, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &a);

  Mat aTa;
  MatTransposeMatMult(a, a, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &aTa);


  Vec d;
  VecDuplicate(dtau, &d);
  VecScale(subDist, 1);
  MatMultTranspose(a, subDist, d);

  Vec dd;
  Vec dd2;
  VecDuplicate(dtau, &dd);
  VecDuplicate(dtau, &dd2);

  if (l1_ != 0) {
    VecSet(dd, l1_);
    MatDiagonalSet(aTa, dd, ADD_ALL_VALUES);
  }

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
  KSPSolve(ksp, d, dtau);

  // -- -- -- -- -- -- -- -- -- --- -- --- ---

  PetscScalar *t;
  PetscScalar *t2;

  VecGetArray(dtau, &t);
  VecGetArray(ti_, &t2);

  for (int i = 0; i < numElementOfInterestIndices; i++) {
    TFloat val = activeStress_->Get(lastTime_, elementsOfInterestMapping_[i]);

    if (t[i] > 10000)
      t[i] = 10000;
    else if (t[i] < -10000)
      t[i] = -10000;

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

  VecRestoreArray(dtau, &t);
  VecRestoreArray(ti_, &t2);

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
  for (auto &i : solidElements_) {
    if (dynamic_cast<CBElementSolidT10T4 *>(i) != 0)
      dynamic_cast<CBElementSolidT10T4 *>(i)->SwitchToT10();
    if (dynamic_cast<CBElementSolidT10RIT4 *>(i) != 0)
      dynamic_cast<CBElementSolidT10RIT4 *>(i)->SwitchToT10();
  }
  lastEstimatorTimeStep_ = time;
  lastTime_ = time;
  return CBStatus::SUCCESS;
} // CBSolverActiveStressEstimator::EstimatorStep

void CBSolverActiveStressEstimator::GenerateElementLaplacian() {
  Mat cMat, subCMat;
  Vec diag;
  Vec b;

  Petsc::CreateSeqMatrix(GetNumberOfSolidElements(), GetNumberOfSolidElements(), 100, &cMat);
  Petsc::CreateSeqVector(GetNumberOfSolidElements(), &diag);
  VecSet(diag, 0);

  for (auto &e : solidElements_) {
    std::unordered_set<TInt> neighbors = adapter_->GetSolver()->GetModel()->GetSolidElementNeighborsCommonFace(
      e->GetIndex());
    for (auto &n : neighbors)
      MatSetValue(cMat, e->GetIndex(), n, -1, INSERT_ALL_VALUES);
  }
  MatDiagonalSet(cMat, diag, INSERT_ALL_VALUES);
  VecDuplicate(diag, &b);
  VecSet(diag, 1);
  MatMult(cMat, diag, b);
  VecScale(b, -1);

  MatDiagonalSet(cMat, b, INSERT_ALL_VALUES);

  MatAssemblyBegin(cMat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(cMat, MAT_FINAL_ASSEMBLY);

  MatGetSubMatrix(cMat, elementsOfInterestIndices_, elementsOfInterestIndices_, MAT_INITIAL_MATRIX, &subCMat);
  MatCreateSeqAIJ(Petsc::Comm(), numElementOfInterestIndices, numElementOfInterestIndices, 0, 0, &elementLaplacian_);
  MatDuplicate(subCMat, MAT_DO_NOT_COPY_VALUES, &elementLaplacian_);
  MatCopy(subCMat, elementLaplacian_, DIFFERENT_NONZERO_PATTERN);
  MatTransposeMatMult(elementLaplacian_, elementLaplacian_, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &lTl_);
  MatDestroy(&cMat);
  MatDestroy(&subCMat);

  VecDestroy(&diag);
  VecDestroy(&b);
} // CBSolverActiveStressEstimator::GenerateElementLaplacian
