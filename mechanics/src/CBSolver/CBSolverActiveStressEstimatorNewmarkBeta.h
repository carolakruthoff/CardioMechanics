//
//  CBSolverActiveStressEstimatorNewmarkBeta.h
//  CardioMechanics_unstable
//
//  Created by Thomas Fritz on 07.06.12.
//  Copyright (c) 2012 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
//

#ifndef CB_SOLVER_ACTIVE_STRESS_ESTIMATOR_NEWMARK_BETA_H
#define CB_SOLVER_ACTIVE_STRESS_ESTIMATOR_NEWMARK_BETA_H

#include <cstdio>

#include "CBSolver.h"
#include "CBSolverEquilibrium.h"
#include "CBSolverNewmarkBeta.h"
#include "CBContactHandling.h"
#include "CBDetermineNodalForces.h"

class CBSolverActiveStressEstimatorNewmarkBeta : public CBSolverNewmarkBeta


{
public:
    
    CBSolverActiveStressEstimatorNewmarkBeta() : CBSolverNewmarkBeta() {}
    std::string GetType(){return("Active Stress Estimator (Newmark Beta)"); }
    CBStatus SolverStep(PetscScalar time, bool forceJacobianAndDampingRecalculation = false);
    virtual void Init(ParameterMap* parameters, CBModel* model);
protected:
private:
    CBStatus EstimatorStep(PetscScalar time);
    CBDataCtrl* activeStress_;
    CBContactHandling* contact_;
    CBDetermineNodalForces* nForces_;
    void GenerateElementLaplacian();
    double currentTime_ = 0;
    double lastTime_=0;
    IS masterNodesIndices_;
    IS nodesOfInterestIndices_;
    IS elementsOfInterestIndices_;
    TInt numMasterNodesIndices_=0;
    TInt numNodesOfInterestIndices_=0;
    TInt numElementOfInterestIndices=0;
    TInt* elementsOfInterestMapping_;
    TInt* masterNodesIndicesNodesOfInterestMapping_;
    TInt* masterNodesIndicesMapping_;
    std::vector<TInt> mat_;
    Vec tau_;
    Vec dx_;
    Vec dist_;
    Vec tmpNodes_;
    Mat elementLaplacian_;
    Vec lastSolution_;
    Vec tmpSolution_;
    TFloat l1_;
    TFloat l2_;
    TFloat l3_;
    int thikonov_ = 1;
    bool staticEstimation_=true;
};

#endif



