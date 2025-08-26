/*
 *  CBRegistrationICS.h
 *  CardioMechanics
 *
//  Created by Thomas Fritz and Emily Reid on 06.09.11. 
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#ifndef CB_REGISTRATION_ITERATIVE_CLOSEST_SURFACE
#define CB_REGISTRATION_ITERATIVE_CLOSEST_SURFACE

#include <map>
#include <vector>
#include "Matrix3.h"

#include <petscksp.h>

#include "CBSolverPlugin.h"
#include "CBCavity.h"
#include "CBElementCavity.h"
#include "CBRegistration.h"

using namespace math_pack;

class CBRegistrationICS : public CBRegistration
{
public:
    CBRegistrationICS() :
        alpha_(0),
        beta_(0),
        minDist_(0),
        status_(CBStatus::REPEAT)
    {}

    virtual ~CBRegistrationICS(){}

    void Init();
    void Apply(TFloat time);
    void ApplyToNodalForces();
    void ApplyToNodalForcesJacobian();
    void FixSourceNodes();

    TFloat GetAverageRegistrationForce(){return(avgRegForce_);}
    CBStatus GetStatus(){return(status_);}
    std::string GetName(){return("NonRigidICS");}

protected:
private:
    void LoadTargetSurfaces(const std::string& nodesFilename, const std::string& surfacesFilename);
    void WriteSurfaceToFile(const std::string& filename);
    void LoadSourceNodesIndices(const std::string& filename);

    std::map<TInt, std::vector<Triangle<TFloat> > > targetSurfaces_;
    std::map<TInt, std::vector<Triangle<TFloat> > > closestTriangle_;
    std::map<TInt, std::vector<bool> >              canBeProjected_;
    std::map<TInt, std::vector<TInt> >              sourceNodesCoordsLocalIndices_;

    TFloat           alpha_;
    TFloat           beta_;
    TFloat           minDist_;
    CBStatus status_;
    TFloat           atol_;
    bool*            fix_;
    TFloat           avgRegForce_;
    typedef CBSolverPlugin   Base;
};
#endif
