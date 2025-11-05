/**@file CBSolverActiveStressTensorEstimator.cpp
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

//  CBSolverActiveStressTensorEstimator.cpp
//  CardioMechanics_unstable
//
//  Created by Thomas Fritz on 07.06.12.
//  Copyright (c) 2012 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
//

#include <stdlib.h>
#include <algorithm>
#include <iterator>



#include "petscversion.h"
#ifndef MY_PETSC_VERSION
# define MY_PETSC_VERSION (PETSC_VERSION_MAJOR*10000 + PETSC_VERSION_MINOR*100 + PETSC_VERSION_SUBMINOR)
#endif // ifndef MY_PETSC_VERSION
#if (MY_PETSC_VERSION >= 30600)
# include <petsc/private/matimpl.h> // for MatAXPY (sparse version of MatAXPY) -- petsc >= 3.6
#else // if (MY_PETSC_VERSION >= 30600)
# include <petsc-private/matimpl.h> // for MatAXPY (sparse version of MatAXPY) -- petsc <= 3.5
#endif // if (MY_PETSC_VERSION >= 30600)

#include "CBSolverActiveStressTensorEstimator.h"
#include "CBContactHandling.h"
#include "CBElementSolidT10T4.h"
#include "CBElementSolidT10RIT4.h"

/// initializes the inverse solver
void CBSolverActiveStressTensorEstimator::Init(ParameterMap *parameters, CBModel *model) {
  // call init function from parent class
  CBSolverEquilibrium::Init(parameters, model);

  // initialize active stress object
  activeStress_ = new CBDataCtrl;
  activeStress_->Init(GetNumberOfElements());
  SetActiveStressDataSource(activeStress_);

  // allocate vectors for displacement
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
  // eki: Wo liegt der Unterschied? elements und threeElements? Eines für nur-stress, das andere fuer die stress+Fasern. Ist der erste damit ueberfluessig?
  std::set<TInt> elementsOfInterest;
  std::set<TInt> threeElementsOfInterest; // diff to CBSolverActiveStressEstimator

  // switch to element type T4 if necessary
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

    // diff to CBSolverActiveStressEstimator (next 5 lines)
    // get indices for threeElements of interest
    auto index = i->GetIndex();
    elementsOfInterest.insert(index);
    threeElementsOfInterest.insert(3*index);
    threeElementsOfInterest.insert(3*index+1);
    threeElementsOfInterest.insert(3*index+2);
  }
  // there are 3 values for each node of interest -> n: amount of evaluation point components
  numNodesOfInterestIndices_ = 3*nodesOfInterest.size();
  numElementOfInterestIndices = elementsOfInterest.size();
  numThreeElementOfInterestIndices = threeElementsOfInterest.size(); // diff to CBSolverActiveStressEstimator

  PetscInt *ni = new PetscInt[numNodesOfInterestIndices_];
  elementsOfInterestMapping_ = new PetscInt[numElementOfInterestIndices];
  // diff to CBSolverActiveStressEstimator (next 6 lines)
  // allocate matrices
  threeElementsOfInterestMapping_ = new PetscInt[numThreeElementOfInterestIndices];
  // maps indices of elements of interest to corresponding solid element indices
  int cc = 0;
  elementOfInterestIndexToSolidElementIndexMapping_ = new TInt[numElementOfInterestIndices];
  DCCtrl::print << "numElementOfInterestIndices " << numElementOfInterestIndices << "\n";
  // diff to CBSolverActiveStressEstimator (next 7 lines)
  for (auto &i : solidElements_) {
    if (elementsOfInterest.find(i->GetIndex()) != elementsOfInterest.end()) {
      elementOfInterestIndexToSolidElementIndexMapping_[cc] = i->GetIndex();
      cc++;
    }
  }
  DCCtrl::print << "cc " << cc << "\n"; // probably for debugging

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

  // thikonov_ = parameters_->Get<TInt>("Solver.ActiveStressEstimator.Thikonov",0);
  // if(thikonov_ != 0 && thikonov_ != 1 && thikonov_ != 2 )
  //     throw std::runtime_error("void CBSolverActiveStressTensorEstimatorNewmarkBeta::Init(ParameterMap* parameters, CBModel* model");

  // li_ are penalty parameters for the thikonov regularization, they can all be set to 0
  // it is advised to set 1 to 2 of these parameters unequal to 0, but never all
  // these parameters unequal to zero are defined in the .xml script
  l1_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l1", 0);
  l2_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l2", 0);
  l3_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l3", 0);
  l4_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l4", 0);
  l5_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l5", 0);
  l6_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l6", 0);
  l7_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l7", 0);  // currently unused in tensor estimator

  // diff to CBSolverActiveStressEstimator (next 11 lines)
  // more parameters, probably for fiber regularisation, currently unused
  l8_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l8", 0);
  l9_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l9", 0);
  l10_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l10", 0);
  l11_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l11", 0);
  l12_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l12", 0); // penalty for large fiber rotation velocity (absolute)
  l13_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l13", 0);
  l14_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l14", 0);
  l15_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l15", 0);
  l16_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l16", 0); // laplace regularization of fiber rotation velocity
  l17_ = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.l17", 0);

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

  // diff to CBSolverActiveStressEstimator (next 6 lines)
  // map local to global indices for three elements of Interest
  it = threeElementsOfInterest.begin();
  for (int i = 0; i < threeElementsOfInterest.size(); i++) {
    // eki: Was sind threeElementsOfInterest ? (fuer Submatrix) Jetzt drei Werte statt einem: stress + phi + theta
    threeElementsOfInterestMapping_[i] = *it;
    it++;
  }

  // Create PETSc index sets (IS) for master nodes, elements and nodes of interest (efficient data access)
  ISCreateGeneral(
    Petsc::Comm(), numMasterNodesIndices_, masterNodesIndicesMapping_, PETSC_COPY_VALUES, &masterNodesIndices_);
  ISCreateGeneral(
    Petsc::Comm(), numElementOfInterestIndices, elementsOfInterestMapping_, PETSC_COPY_VALUES,
    &elementsOfInterestIndices_);
  // diff to CBSolverActiveStressEstimator (next 3 lines)
  ISCreateGeneral(
    Petsc::Comm(), numThreeElementOfInterestIndices, threeElementsOfInterestMapping_, PETSC_COPY_VALUES,
    &threeElementsOfInterestIndices_);

  ISCreateGeneral(Petsc::Comm(), numNodesOfInterestIndices_, ni, PETSC_COPY_VALUES, &nodesOfInterestIndices_);
  // duplicates nodes_ into dist_ -> storing distances or residuals
  VecDuplicate(nodes_, &dist_);
  int counter = 0;   // diff to CBSolverActiveStressEstimator
  // Switch to T10 elements, this is usually done at the end of the estimator step
  // since the estimation is done with T4 elements, independent of the type of mesh elements
  for (auto &i : solidElements_) {
    counter++;   // diff to CBSolverActiveStressEstimator
    if (dynamic_cast<CBElementSolidT10T4 *>(i) != 0)
      dynamic_cast<CBElementSolidT10T4 *>(i)->SwitchToT10();
    if (dynamic_cast<CBElementSolidT10RIT4 *>(i) != 0)
      dynamic_cast<CBElementSolidT10RIT4 *>(i)->SwitchToT10();
  }
  DCCtrl::print << "switched " << counter << " elements to T4\n";

  // set up laplacian matrix (spatial smoothing regularization)
  GenerateElementLaplacian();
  // creates and duplicates Petsc vectors
  // TODO: Ist ti_ die Loesung aus dem letzten Zeitschritt? -> fuer zeitliche Regularisierung? aktueller Zeitschritt, ti_ -> T_{i}
  GenerateElementLaplacianFibers();
  Petsc::CreateSeqVector(numThreeElementOfInterestIndices, &ti_);  // diff to CBSolverActiveStressEstimator
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
void CBSolverActiveStressTensorEstimator::UpdateMasterNodesOfInterest() {
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
} // CBSolverActiveStressTensorEstimator::UpdateMasterNodesOfInterest

/// eki: the actual optimization / parameter estimation step, calls estimator + solver + processes the results
CBStatus CBSolverActiveStressTensorEstimator::SolverStep(PetscScalar time) {
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
  // in ActiveStressEstimator alpha is read at this point, but never used
  //    double beta = parameters_->Get<TFloat>("Solver.ActiveStressEstimator.Beta",1e8);


  // PetscScalar* t;
  // VecGetArray(lastSolution_, &t);

  //    for(int i=0; i < numElementOfInterestIndices ; i++)
  //    {
  //        activeStress_->Set(time,elementsOfInterestMapping_[i],t[i] - dTau_[i]);
  //        dTau_[i] = 0;
  //    }

  // VecRestoreArray(lastSolution_, &t);

  //    contact_->SetAlpha(alpha);

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

    ////////////
    // Estimate: compute as + fiber angles
    ///////////

    // do estimator step
    EstimatorStep(time, i);

    //        contact_->SetAlpha(contact_->GetAlpha()/100);
    contact_->Apply(time); // apply from ContactHandling
    UpdateGhostNodesAndLinkToAdapter();             // Was sind die Ghost nodes? Was machen sie? -> Kommunikation am Rand (Parallelisierung)
    CreateNodesJacobianAndLinkToAdapter();
    UpdateActiveStress(time);

    /////////////////
    // SolverStep: apply estimated values and test if solution is good
    /////////////////

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
  if (rc == CBStatus::FAILED) { // Schrittweite zu gross?
    PetscScalar *t;
    VecGetArray(ti1_, &t);

    for (int i = 0; i < numElementOfInterestIndices; i++) {
      // TODO: Werden die Ableitungen berechnet?
      // TODO: Zshg dTau, lastDtau <-> ti, ti1, ti2 (T_i, T_{i-1}, T_{i-2}) ?
      // TODO: Wo werden die Werte fuer dTau oder dt berechnet?
      //
      activeStress_->Set(time, elementsOfInterestMapping_[i], -t[3*i]); // 3*i: diff to CBSolverActiveStressEstimator

      // Fasern wieder zurueckdrehen ( in t[3*i+1] und t[3*i+2] ), oder? Stand bereits da :-)
      // TODO: Funktioniert das reversibel? -> Wahrscheinlich nicht -> in zwei Schritten machen
      solidElements_[elementOfInterestIndexToSolidElementIndexMapping_[i]]->RotateFiberBy(0, t[3*i+2]); // diff to CBSolverActiveStressEstimator
      solidElements_[elementOfInterestIndexToSolidElementIndexMapping_[i]]->RotateFiberBy(t[3*i+1], 0); // diff to CBSolverActiveStressEstimator
    }
    VecRestoreArray(ti1_, &t); // pointer t is no langer valid after this call, I'm done accessing ti1_
    VecCopy(ti1_, ti_);
  } else {   // Schaetzer hat funktioniert
    VecCopy(ti1_, ti2_);

    // fetch array from ti1_ as t
    PetscScalar *t;
    VecGetArray(ti1_, &t);

    for (int i = 0; i < numElementOfInterestIndices; i++) {
      // TODO: !!! Hier muss noch etwas mit t[3*i+1] und t[3*i+2] geschehen?
      t[3*i] = -activeStress_->Get(time, elementsOfInterestMapping_[i]); // 3*i: diff to CBSolverActiveStressEstimator
    }
    VecRestoreArray(ti1_, &t);
    VecCopy(ti1_, ti_);

    // diff to CBSolverActiveStressEstimator (next 7 lines)
    // Winkelaenderung mit 0 initialisieren -> hat irgendwie nicht funktioniert (?)
    VecGetArray(ti_, &t);
    for (int i = 0; i < numElementOfInterestIndices; i++) {
      t[3*i+1] = 0;
      t[3*i+2] = 0;
    }
    VecRestoreArray(ti_, &t);
  }

  return rc;
} // CBSolverActiveStressTensorEstimator::SolverStep

CBStatus CBSolverActiveStressTensorEstimator::EstimatorStep(PetscScalar time, int step) {
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
  Petsc::CreateSeqVector(numThreeElementOfInterestIndices, &dtau); // 'numThreeElement': diff to CBSolverActiveStressEstimator
  Petsc::CreateSeqMatrix(3*numNodes_, 3*GetNumberOfElements(), 5000, &dfdtau); // 3*i: diff to CBSolverActiveStressEstimator

  ////
  // determine parameters, compute the derivatives for active stress and fiber rotation
  ////

  adapter_->LinkNodalForcesActiveStressTensorAndFiberOrientationJacobian(dfdtau); // diff to CBSolverActiveStressEstimator
  formulation_->CalcNodalForcesActiveStressTensorAndFiberOrientationJacobian(); // diff to CBSolverActiveStressEstimator
  //formulation is an object ob CBSolverActiveStressEstimator
  // stores tangential stiffness matrix in dfdtau
  MatAssemblyBegin(dfdtau, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(dfdtau, MAT_FINAL_ASSEMBLY);
  // extract submatrix of tangential stiffness matrix containing only nodes of interest
  MatGetSubMatrix(dfdtau, nodesOfInterestIndices_, threeElementsOfInterestIndices_, MAT_INITIAL_MATRIX, &subdfdtau); // 'ThreeElement': diff to CBSolverActiveStressEstimator

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

    //        DCCtrl::print << "f >>>\n";
    //        PetscViewer viewer;
    //        PetscViewerASCIIOpen(PETSC_COMM_WORLD, "./matrix.txt", &viewer);
    //        PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_DENSE);
    //        MatView(f, viewer);
    //        PetscViewerDestroy(&viewer);
    //        DCCtrl::print << "<<< f\n";

    MatLUFactorSymbolic(f, subdfdx, perm, iperm, 0);               // ### LU Zerlegung von subdfdx vorbereiten, Ergebnis in f
    MatLUFactorNumeric(f, subdfdx, 0);                               // LU-Zerlegung, invertieren

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
      VecSetValue(ba, masterNodesIndicesNodesOfInterestMapping_[i], 1, INSERT_VALUES);  // ### rechte Seite um reduzierte Inverse zu bestimmen

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

    // check if inv_ is a nullpointer -> if yes create a matrix inv_ as duplicate of inv
    // but do not copy the values from inv
    if (inv_ == 0)
      MatDuplicate(inv, MAT_DO_NOT_COPY_VALUES, &inv_);

    // assemble inverse matrix in inv_
    MatCopy(inv, inv_, DIFFERENT_NONZERO_PATTERN);
    MatAssemblyBegin(inv_, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(inv_, MAT_FINAL_ASSEMBLY);
    // reuse previous inverse if applicable
  } else {
    Petsc::debug << "Reusing previous matrix\n";
    MatCopy(inv_, inv, DIFFERENT_NONZERO_PATTERN);
    MatAssemblyBegin(inv, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(inv, MAT_FINAL_ASSEMBLY);
  }

  // compute a = inv * subdfdtau = L * k^_d
  MatMatMult(inv, subdfdtau, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &a);
  //compute aTa = a^T * a

  // hier ist alles dreimal so gross, wg Definition von subdfdtau
  Mat aTa;
  DCCtrl::print << "Compute aTa...\n";
  clock_t startTime = clock();
  MatTransposeMatMult(a, a, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &aTa);  // ### A^T*A, bzw. die linke Seite um Ableitung=0 zu bestimmen
  DCCtrl::print << "Compute aTa... OK, took " << double(clock() - startTime) / (double)CLOCKS_PER_SEC<< " seconds.\n";

  // Apply tikhonov regularization terms to stabilize solution
  // each regularization term adds a weighted identity or Laplacian matrix to aTa and adjusts d
  Vec d; // right-hand side vector
  VecDuplicate(dtau, &d);
  VecScale(subDist, 1);
  MatMultTranspose(a, subDist, d);  // ### d = A^T g, bzw. die rechte Seite um Ableitung 0 zu bestimmen

  Vec dd; // temporary vector
  // in ActiveStressEstimator there is also a dd2
  VecDuplicate(dtau, &dd);

  // if(thikonov_ == 0)
  // {
  //     VecSet(dd, l1_);
  //     MatDiagonalSet(aTa, dd, ADD_ALL_VALUES);  // aTa = aTa + dd
  // }
  // else if(thikonov_ == 1 || thikonov_ == 2)
  // {
  //     VecSet(dd, l2_*Base::GetTimeStepSize());
  //     MatDiagonalSet(aTa, dd, ADD_ALL_VALUES);           // aTa = aTa + dd
  //     MatAXPY(aTa, l3_, lTl_,DIFFERENT_NONZERO_PATTERN); // aTa = l3*lTl + aTa
  // }
  // else
  //     throw std::runtime_error("CBStatus CBSolverActiveStressTensorEstimatorNewmarkBeta::EstimatorStep(PetscScalar time): Thikonov must be 0 [L2-Norm],1 [L2-Norm Laplace dtau] or 2 [L2-Norm Laplace tau]\n");

  /*
     // alte Version der Tikhonov-Regularisierung
     Vec t1;
     Vec t2;
     VecDuplicate(d, &t1);
     VecDuplicate(d, &t2);
     PetscScalar* v;

     VecGetArray(t1,&v);

     for(int i=0; i<numThreeElementOfInterestIndices;i++)
     v[i] = dTau_[i];

     VecRestoreArray(t1, &v);

     // if(thikonov_ == 2);
     //     VecAXPY(t1, 1, lastSolution_);

     MatMult(lTl_,t1,t2);
     VecScale(t2,l3_);

     VecAXPY(d, -1.0, t2);
     VecAXPY(d, l2_, lastDtau_);

     VecDestroy(&t1);
     VecDestroy(&t2);
   */


  // diff to ActiveStressEstimator
  // useful indices for setting values during regularization
  PetscInt *stressIndices;
  PetscInt *phiIndices;
  PetscInt *thetaIndices;
  stressIndices = new PetscInt[numElementOfInterestIndices];
  phiIndices = new PetscInt[numElementOfInterestIndices];
  thetaIndices = new PetscInt[numElementOfInterestIndices];
  for (int i = 0; i < numElementOfInterestIndices; i++) {
    stressIndices[i] = 3*i+0;
    phiIndices[i]    = 3*i+1;
    thetaIndices[i]  = 3*i+2;
  }


  // TODO: Das leider 1 zu 1 die Regularisierungsterme fuer die aktiven Kraefte, aber Fasern brauchen evtl etwas anderes
  // -> nur jeden dritten Eintrag setzen?
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

  // diff to ActiveStressEstimator (if for l12 and l16 as well as the 3 deletes afterwards)
  // regularization for minimal d-phi, d-theta (minimal rotation change)
  if (l12_ != 0) {
    VecSet(dd, 0);
    PetscScalar *values = new PetscScalar[numElementOfInterestIndices];
    for (int i = 0; i < numElementOfInterestIndices; i++) {
      values[i] = l12_;
    }
    VecSetValues(dd, numElementOfInterestIndices, phiIndices, values, INSERT_VALUES);
    VecSetValues(dd, numElementOfInterestIndices, thetaIndices, values, INSERT_VALUES);

    MatDiagonalSet(aTa, dd, ADD_ALL_VALUES);
    VecAXPY(d, -l12_, ti_);
    delete[] values;
  }

  // regularization with laplace for d-phi, d-theta (difference, not exact value!)
  if (l16_ != 0) {
    MatAXPY(aTa, l16_, lTlFibers_, DIFFERENT_NONZERO_PATTERN);
    MatMult(lTlFibers_, ti_, dd);
    VecAXPY(d, -l16_, dd);

    // MatView(aTa, PETSC_VIEWER_STDOUT_WORLD);
    // MatView(lTl_, PETSC_VIEWER_STDOUT_WORLD);
  }

  delete stressIndices;
  delete phiIndices;
  delete thetaIndices;


  // -- -- -- -- -- -- -- -- -- --- -- --- ---
  DCCtrl::print << "ksp...\n";
  startTime = clock();
  // Set up linear solver (KSP) with direct LU factorization
  // uses MUMPS
  KSP ksp;
  KSPCreate(Petsc::Comm(), &ksp);
  PC pc;
  KSPGetPC(ksp, &pc);
  KSPSetType(ksp, "preonly");
  PCFactorSetMatSolverPackage(pc, "mumps");
  KSPSetFromOptions(ksp);

  // PCSetType(pc, PCLU);
  PCSetType(pc, PCCHOLESKY); // diff type than in ActiveStressEstimator
  // Cholesky factorization is specialized for symmetric positive definite matrices

  KSPSetOperators(ksp, aTa, aTa);

  // TODO: !!! Das dauert lange, koennte man das vlt nur ein einziges Mal initialisieren und spaeter die Matrizen tauschen?
  DCCtrl::print << "ksp setup...\n";
  KSPSetUp(ksp);
  DCCtrl::print << "ksp solve...\n";
  // solves linear system aTa * dtau = d to find changes in active stress (dtau)
  KSPSolve(ksp, d, dtau);
  DCCtrl::print << "ksp... OK, took " << double(clock() - startTime) / (double)CLOCKS_PER_SEC<< " seconds.\n";

  // -- -- -- -- -- -- -- -- -- --- -- --- ---

  PetscScalar *t;
  PetscScalar *t2;

  VecGetArray(dtau, &t);
  VecGetArray(ti_, &t2);
  // update active stress values for each element of interest (potentially have to edit this?!)
  // Beschraenkung der maximalen AS-Werte auf 10.000
  for (int i = 0; i < numElementOfInterestIndices; i++) {
    TFloat val = activeStress_->Get(lastTime_, elementsOfInterestMapping_[i]);

    // diff to ActiveStressEstimator: all factors 3*i for the following if conditions
    if (t[3*i] > 10000)
      t[3*i] = 10000;
    else if (t[3*i] < -10000)
      t[3*i] = -10000;
    // applies bounds to change in stress to ensure numerical stability
    if ((val - t[3*i]) < 0) {
      t[3*i] = val;
      activeStress_->Set(time, elementsOfInterestMapping_[i], val - t[3*i]);
    } else if ((val - t[3*i]) > 2e5) {
      t[3*i] = val - 2e5;
      activeStress_->Set(time, elementsOfInterestMapping_[i], val - t[3*i]);
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
      activeStress_->Set(time, elementsOfInterestMapping_[i], val - t[3*i]); // 3*i: diff to ActiveStressEstimator
    }

    // TODO: dTau_ oder t2?
    // TODO: !!! wird letztlich mit t (dtau) oder t2 (ti) gerechnet?
    // TODO: muss hier dTau_ noch ersetzt werden?
    t2[3*i] += t[3*i]; // 3*i: diff to ActiveStressEstimator

    // dTau_[3*i] += t[3*i];
    // dTau_[3*i+1] += t[3*i+1];
    // dTau_[3*i+2] += t[3*i+2];

    t2[3*i+1] = t[3*i+1]; // diff to ActiveStressEstimator
    t2[3*i+2] = t[3*i+2]; // diff to ActiveStressEstimator

    solidElements_[elementOfInterestIndexToSolidElementIndexMapping_[i]]->RotateFiberBy(-t[3*i+1], 0);
    solidElements_[elementOfInterestIndexToSolidElementIndexMapping_[i]]->RotateFiberBy(0, -t[3*i+2]);
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
} // CBSolverActiveStressTensorEstimator::EstimatorStep

/// eki: fills elementLaplacian_ with entries, actually just a connectivity matrix for regularization
void CBSolverActiveStressTensorEstimator::GenerateElementLaplacian() {
  Mat cMat, subCMat; // C
  Vec diag;
  Vec b;
  // create sequential sparse matrix (cMat) with dimensions equal to number of solid elements
  // 3*: diff to ActiveStressEstimator (in next 2 lines)
  Petsc::CreateSeqMatrix(3*GetNumberOfSolidElements(), 3*GetNumberOfSolidElements(), 100, &cMat);  // C-Matrix, connectivity
  Petsc::CreateSeqVector(3*GetNumberOfSolidElements(), &diag);                                  // diag-Matrix
  // create vector of same size as cMat, initialize to 0
  VecSet(diag, 0);

  int max = 0; // diff to ActiveStressEstimator
  // iterate over all solid elements
  // build connectivity matrix where each element is connected to its neighbors
  for (auto &e : solidElements_) {
    int counter = 0; // diff to ActiveStressEstimator
    // retrieve neighbor for each element (share common face)
    std::unordered_set<TInt> neighbors = adapter_->GetSolver()->GetModel()->GetSolidElementNeighborsCommonFace(
      e->GetIndex());
    // sets value of cMat at position (e->GetIndex(), n) to -1 for each neighbor n
    for (auto &n : neighbors) {
      // 3*: diff to ActiveStressEstimator (in next line)
      MatSetValue(cMat, 3*e->GetIndex(), 3*n, -1, INSERT_ALL_VALUES);                        // C-Matrix
      // MatSetValue(cMat,3*(*e)->GetIndex()+1,3*(*n)+1,-1,INSERT_ALL_VALUES);
      // MatSetValue(cMat,3*(*e)->GetIndex()+2,3*(*n)+2,-1,INSERT_ALL_VALUES);
      counter++; // diff to ActiveStressEstimator
    }

    // DCCtrl::print << "numNeighbors: " << counter << "\n";
    if (counter > max) max = counter; // diff to ActiveStressEstimator
  }
  DCCtrl::print << "max numNeighbors: " << max << "\n"; // diff to ActiveStressEstimator

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
  MatGetSubMatrix(cMat, threeElementsOfInterestIndices_, threeElementsOfInterestIndices_, MAT_INITIAL_MATRIX, &subCMat); // 'threeElements': diff to ActiveStressEstimator

  // next commented part is used in ActiveStressEstimator
  // MatCreateSeqDense(Petsc::Comm(), numThreeElementOfInterestIndices, numThreeElementOfInterestIndices, 0, &elementLaplacian_);
  // ------------------
  // MatCreateSeqAIJ(Petsc::Comm(), numThreeElementOfInterestIndices, numThreeElementOfInterestIndices, 0, 0, &elementLaplacian_);
  // MatDuplicate(subCMat, MAT_DO_NOT_COPY_VALUES, &elementLaplacian_);
  // ------------------
  // MatCopy(subCMat, elementLaplacian_, DIFFERENT_NONZERO_PATTERN);

  // MatDuplicate(subCMat, MAT_DO_NOT_COPY_VALUES, &lTl_);

  DCCtrl::print << "############################### generate elem Laplacian: \n"; // diff to ActiveStressEstimator

  // MatTransposeMatMult(elementLaplacian_, elementLaplacian_,MAT_INITIAL_MATRIX , PETSC_DEFAULT, &lTl_);
  // MatView(lTl_, PETSC_VIEWER_STDOUT_WORLD);


  Mat spLaplacian;
  // create and duplicate element laplacian (copy subCMAt into it)
  MatCreateSeqAIJ(
    Petsc::Comm(), numThreeElementOfInterestIndices, numThreeElementOfInterestIndices, 20, NULL, &spLaplacian); // diff to ActiveStressEstimator
  MatCopy(subCMat, spLaplacian, DIFFERENT_NONZERO_PATTERN); // diff to ActiveStressEstimator
  Mat splTl;  // diff to ActiveStressEstimator
  MatCreateSeqAIJ(Petsc::Comm(), numThreeElementOfInterestIndices, numThreeElementOfInterestIndices, 20, NULL, &splTl); // diff to ActiveStressEstimator
  // represents stiffness-like operator)
  MatTransposeMatMult(spLaplacian, spLaplacian, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &splTl); // diff to ActiveStressEstimator

  MatCreateSeqAIJ(Petsc::Comm(), numThreeElementOfInterestIndices, numThreeElementOfInterestIndices, 20, NULL, &lTl_);  // diff to ActiveStressEstimator
  MatCopy(splTl, lTl_, DIFFERENT_NONZERO_PATTERN);  // diff to ActiveStressEstimator


  DCCtrl::print << "############################### generate elem Laplacian: OK\n";  // diff to ActiveStressEstimator
  // cleanup
  MatDestroy(&cMat);
  MatDestroy(&subCMat);
  MatDestroy(&splTl);
  MatDestroy(&spLaplacian);

  VecDestroy(&diag);
  VecDestroy(&b);
} // CBSolverActiveStressTensorEstimator::GenerateElementLaplacian

/// fills elementLaplacianFibers_ with entries, actually just a connectivity matrix for Fiber regularization (not existant in ActiveStressEstimator)
void CBSolverActiveStressTensorEstimator::GenerateElementLaplacianFibers() {
  Mat cMat, subCMat; // C
  Vec diag;
  Vec b;

  Petsc::CreateSeqMatrix(3*GetNumberOfSolidElements(), 3*GetNumberOfSolidElements(), 100, &cMat);  // C-Matrix, connectivity
  Petsc::CreateSeqVector(3*GetNumberOfSolidElements(), &diag);                                  // diag-Matrix
  VecSet(diag, 0);

  int max = 0;
  for (auto &e : solidElements_) {
    int counter = 0;
    std::unordered_set<TInt> neighbors = adapter_->GetSolver()->GetModel()->GetSolidElementNeighborsCommonFace(
      e->GetIndex());
    for (auto &n : neighbors) {
      // MatSetValue(cMat,3*(*e)->GetIndex(),3*(*n),-1,INSERT_ALL_VALUES);                        // C-Matrix for tension
      MatSetValue(cMat,
                  3*e->GetIndex()+1,
                  3*n+1,
                  -1,
                  INSERT_ALL_VALUES);                    // C-Matrix for phi
      MatSetValue(cMat,
                  3*e->GetIndex()+2,
                  3*n+2,
                  -1,
                  INSERT_ALL_VALUES);                    // C-Matrix for theta
      counter++;
    }

    // DCCtrl::print << "numNeighbors: " << counter << "\n";
    if (counter > max) max = counter;
  }
  DCCtrl::print << "max numNeighbors: " << max << "\n";

  MatDiagonalSet(cMat, diag, INSERT_ALL_VALUES);               // C = C + diag
  VecDuplicate(diag, &b);              // b gets same size etc as diag, just allocates memory
  VecSet(diag, 1);                     // change all entries to 1
  MatMult(cMat, diag, b);              // b = C * diag
  VecScale(b, -1);                      // b = -b

  MatDiagonalSet(cMat, b, INSERT_ALL_VALUES);  // C = C + b

  MatAssemblyBegin(cMat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(cMat, MAT_FINAL_ASSEMBLY);

  MatGetSubMatrix(cMat, threeElementsOfInterestIndices_, threeElementsOfInterestIndices_, MAT_INITIAL_MATRIX, &subCMat);

  DCCtrl::print << "############################### generate elem Laplacian Fibers: \n";

  Mat spLaplacian;
  MatCreateSeqAIJ(
    Petsc::Comm(), numThreeElementOfInterestIndices, numThreeElementOfInterestIndices, 20, NULL, &spLaplacian);
  MatCopy(subCMat, spLaplacian, DIFFERENT_NONZERO_PATTERN);
  Mat splTl;
  MatCreateSeqAIJ(Petsc::Comm(), numThreeElementOfInterestIndices, numThreeElementOfInterestIndices, 20, NULL, &splTl);

  MatTransposeMatMult(spLaplacian, spLaplacian, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &splTl);


  //    MatView(spLaplacian, PETSC_VIEWER_STDOUT_WORLD);
  //    MatView(splTl, PETSC_VIEWER_STDOUT_WORLD);
  //    MatView(lTl_, PETSC_VIEWER_STDOUT_WORLD);
  MatCreateSeqAIJ(
    Petsc::Comm(), numThreeElementOfInterestIndices, numThreeElementOfInterestIndices, 20, NULL, &lTlFibers_);
  MatCopy(splTl, lTlFibers_, DIFFERENT_NONZERO_PATTERN);

  DCCtrl::print << "############################### generate elem Laplacian Fibers: OK\n";

  MatDestroy(&cMat);
  MatDestroy(&subCMat);
  MatDestroy(&splTl);
  MatDestroy(&spLaplacian);

  VecDestroy(&diag);
  VecDestroy(&b);
} // CBSolverActiveStressTensorEstimator::GenerateElementLaplacianFibers


void CBSolverActiveStressTensorEstimator::UpdateActiveStress(PetscScalar time) {
  // PetscInt cnt = 0;
  // if(activeStressData_) {
  //    for(auto& it : solidElements_) {
  //
  // // todo: this should be able to depend dynamically on everything instead of being set just once, similar to the constitutive model
  //
  //     activeStressData_->SetTime(time);
  //     it->GetMaterial()->GetConstitutiveModel()->SetTemplateForce(activeStressData_->Get(it->GetIndex(), it->GetMaterialIndex()));
  //     // active stress currently gets scaled by the constitutive model, but should rather be in ActiveStressModel
  //
  //     // todo: instead of setting this here to zero, it should be remove from CalcNodalForcesHelperFunction in CBElementT4, CBElementT10, (partially done) and so on...
  //     activeStressTensorComponents_[cnt]        =  activeStressData_->Get(time, it->GetIndex(), it->GetMaterialIndex() );
  //     activeStressTensorComponentsIndices_[cnt] = it->GetLocalIndex() + activeStressLowerIndex_;
  //     cnt++;
  //        }
  //
  // VecSetValues(activeStress_, cnt, activeStressTensorComponentsIndices_, activeStressTensorComponents_, INSERT_VALUES);
  // VecAssemblyBegin(activeStress_);
  // VecAssemblyEnd(activeStress_);
  // }
}
