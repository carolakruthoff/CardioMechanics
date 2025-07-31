//
//  CBSolverActiveStressTensorEstimator.h
//  CardioMechanics_unstable
//
//  Created by Thomas Fritz on 07.06.12.
//  Copyright (c) 2012 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
//

#pragma once

#include <cstdio>

#include "CBSolver.h"
#include "CBSolverEquilibrium.h"
#include "CBSolverNewmarkBeta.h"
#include "CBContactHandling.h"
#include "CBDetermineNodalForces.h"

class CBSolverActiveStressTensorEstimator : public CBSolverEquilibrium


{
public:
    
    CBSolverActiveStressTensorEstimator() : CBSolverEquilibrium() {}
    std::string GetType() override {return("Active Stress Tensor Estimator (Static)"); }
    CBStatus SolverStep(PetscScalar time);
    virtual void Init(ParameterMap* parameters, CBModel* model) override;
protected:
private:
    void UpdateMasterNodesOfInterest();
    CBStatus EstimatorStep(PetscScalar time,int step);
    CBDataCtrl* activeStress_;
    CBContactHandling* contact_;
    CBDetermineNodalForces* nForces_;
    void GenerateElementLaplacian();
    void GenerateElementLaplacianFibers();
    double lastTime_=0;
    IS masterNodesIndices_;
    IS nodesOfInterestIndices_;
    IS elementsOfInterestIndices_;
    IS threeElementsOfInterestIndices_;
    
    TInt numMasterNodesIndices_=0;
    TInt numNodesOfInterestIndices_=0;
    TInt numElementOfInterestIndices=0;
    TInt numThreeElementOfInterestIndices=0;
    
    TInt* elementsOfInterestMapping_;
    TInt* threeElementsOfInterestMapping_;
    
    TInt* masterNodesIndicesNodesOfInterestMapping_;
    TInt* masterNodesIndicesMapping_;
    TInt* elementOfInterestIndexToSolidElementIndexMapping_;
    std::vector<TInt> mat_;
    Mat lTl_=0;
    Mat lTlFibers_=0;
    Vec dx_;
    Vec dist_;
    Vec tmpNodes_;

    //    Vec lastSolution_;            // -> ti1_
    // TFloat* dTau_ = 0;               // -> ti_
    // Vec lastDtau_;                   // -> ti1_
    // int thikonov_ = 1;
    
    Vec ti_;
    Vec ti1_;
    Vec ti2_;


    TFloat l1_ = 1e-55;
    TFloat l2_ = 1e-25;
    TFloat l3_ = 1e-15;
    TFloat l4_ = 0;
    TFloat l5_ = 0;
    TFloat l6_ = 0;
    TFloat l7_ = 0;

    TFloat l8_ = 0;
    TFloat l9_ = 0;
    TFloat l10_ = 0;
    TFloat l11_ = 0;
    TFloat l12_ = 0;
    TFloat l13_ = 0;
    TFloat l14_ = 0;
    TFloat l15_ = 0;
    TFloat l16_ = 0;
    TFloat l17_ = 0;
  
    PetscInt* nim_ = 0;
    std::set<TInt> mn_;
    
    Mat inv_ = 0;
    TFloat lastEstimatorTimeStep_ = -1;
};
