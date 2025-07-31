//
//  CBDetermineNodalForces.h
//  CardioMechanics_unstable
//
//  Created by Thomas Fritz on 09.03.13.
//
//

#ifndef CB_DETERMINE_NODAL_FORCES
#define CB_DETERMINE_NODAL_FORCES 

#include "CBSolverPlugin.h"
#include "CBDataFromFile.h"

class CBDetermineNodalForces : public CBSolverPlugin
{
public:
    void Apply(TFloat time){}
    void SaveNodes();
    std::string GetName(){return(std::string("DetermineNodalForces"));}
    void Init();
    void ApplyToNodalForces();
    void ApplyToNodalForcesJacobian();
protected:
private:
    Vec originalNodes_=0;
    Vec diag_=0;
    PetscScalar k_ = 1e3;
};
#endif
