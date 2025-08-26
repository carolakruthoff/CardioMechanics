/*
 *  CBRegistrationNonRigidICP.h
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 06.09.11.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#ifndef CB_REGISTRATION_ITERATIVE_CLOSEST_POINT
#define CB_REGISTRATION_ITERATIVE_CLOSEST_POINT

#include <map>
#include <vector>
#include "Matrix3.h"

#include <petscksp.h>

#include "CBSolverPlugin.h"
#include "CBSolverPlugin.h"
#include "CBCavity.h"
#include "CBElementCavity.h"
#include "CBRegistration.h"


using namespace math_pack;


class CBRegistrationNonRigidICP : public CBRegistration
{
public:
    CBRegistrationNonRigidICP() : alpha_(0), beta_(0), minDist_(0){}
    virtual ~CBRegistrationNonRigidICP(){}

    void Init();
    void Apply(TFloat time);
    void ApplyToNodalForces();
    void ApplyToNodalForcesJacobian();
    CBStatus GetStatus(){throw std::runtime_error("CBRegistrationNonRigidICP::GetStatus(): This function is not implemented for the requested child class."); return CBStatus::FAILED; }
    std::string GetName(){return "NonRigidICP";}
protected:
private:
    void LoadTargetNodes(const std::string& filename);
    void LoadSourceNodesIndices(const std::string& filename);
    
    std::map<TInt, std::vector<Vector3<TFloat> > > targetNodes_;
    std::map<TInt, std::vector<Vector3<TFloat> > > closestNeighbors_;
    std::map<TInt, std::vector<TInt> >        sourceNodesCoordsLocalIndices_;
    TFloat alpha_;
    TFloat beta_;
    TFloat minDist_;

    typedef CBSolverPlugin   Base;
};


#endif
