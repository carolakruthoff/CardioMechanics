//
//  CBSolverActiveStressEstimator.h
//  CardioMechanics_unstable
//
//  Created by Thomas Fritz on 07.06.12.
//  Copyright (c) 2012 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
//

#ifndef CB_SOLVER_ACTIVE_STRESS_ESTIMATOR_BETA_H
#define CB_SOLVER_ACTIVE_STRESS_ESTIMATOR_BETA_H

#include <cstdio>

#include "CBSolver.h"
#include "CBSolverEquilibrium.h"
#include "CBSolverNewmarkBeta.h"
#include "CBContactHandling.h"
#include "CBDetermineNodalForces.h"

class CBSolverActiveStressEstimator : public CBSolverEquilibrium


{
public:
    
    CBSolverActiveStressEstimator() : CBSolverEquilibrium() {}
    virtual std::string GetType() override {return("Active Stress Estimator (Static)"); }
    // CBStatus SolverStep(PetscScalar time);
    virtual void Init(ParameterMap* parameters, CBModel* model) override;
protected:
    // ck620
    CBStatus SolverStep(PetscScalar time, bool forceJacobianAndDampingRecalculation = false);
private:
    void UpdateMasterNodesOfInterest();
    CBStatus EstimatorStep(PetscScalar time,int step);
    CBDataCtrl* activeStress_;
    CBContactHandling* contact_;
    CBDetermineNodalForces* nForces_;
    void GenerateElementLaplacian();
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
    Mat lTl_=0;
  
    Vec dx_;
    Vec dist_;
    Vec tmpNodes_;
    Mat elementLaplacian_;
    
    Vec ti_;   // T_i
    Vec ti1_;  // T_{i-1}
    Vec ti2_;  // T_{i-2}
  
    TFloat l1_ = 0;
    TFloat l2_ = 0;
    TFloat l3_ = 0;
    TFloat l4_ = 0;
    TFloat l5_ = 0;
    TFloat l6_ = 0;
    TFloat l7_ = 0;
    
    PetscInt* nim_ = 0;
    std::set<TInt> mn_;
    
    Mat inv_ = 0;
    TFloat lastEstimatorTimeStep_ = -1;
};

#endif



