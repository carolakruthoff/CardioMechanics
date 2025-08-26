/*
 *  CBNonRigidICPRegistration.h
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 06.09.11.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#ifndef CB_REGISTRATION
#define CB_REGISTRATION

#include <map>
#include <vector>
#include "Matrix3.h"

#include <petscksp.h>

#include "CBSolverPlugin.h"
#include "CBCavity.h"
#include "CBElementCavity.h"



using namespace math_pack;




class CBRegistration : public CBSolverPlugin
{
public:
    CBRegistration(){}
    virtual ~CBRegistration(){}
    virtual TFloat GetGlobalFitness(){throw std::runtime_error("CBRegistration::GetGlobalFitness(): This function is not implemented for the requested child class."); return -1; }
    virtual void FixSourceNodes(){}
protected:
    virtual void LoadTargetNodes(const std::string& filename){}
    virtual void LoadSourceNodesIndices(const std::string& filename){}
private:
    typedef CBSolverPlugin   Base;
};


#endif
