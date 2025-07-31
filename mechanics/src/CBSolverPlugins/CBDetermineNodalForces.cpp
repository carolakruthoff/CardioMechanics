;//
//  CBDetermineNodalForces.h.cpp
//  CardioMechanics_unstable
//
//  Created by Thomas Fritz on 09.03.13.
//
//

#include <vector>

#include "CBDetermineNodalForces.h"
#include "CBSolver.h"

void CBDetermineNodalForces::Init()
{
    VecDuplicate(adapter_->GetNodes(),&diag_);
    VecSet(diag_,k_);
    VecAssemblyBegin(diag_);
    VecAssemblyEnd(diag_);
}


void CBDetermineNodalForces::SaveNodes()
{
    if(!originalNodes_)
        VecDuplicate(adapter_->GetNodes(),&originalNodes_);
    VecCopy(adapter_->GetNodes(),originalNodes_);
}



void CBDetermineNodalForces::ApplyToNodalForces()
{
    if(!IsActive())
        return;
    Vec f;
    VecDuplicate(originalNodes_,&f);
    VecCopy(originalNodes_,f);
    VecAXPY(f,-1,adapter_->GetNodes());
    VecScale(f,-k_);
    adapter_->AddToForces(f);
    VecDestroy(&f);
}


void CBDetermineNodalForces::ApplyToNodalForcesJacobian()
{
    if(!IsActive())
        return;
    adapter_->AddToNodalForcesJacobianDiagonal(diag_);
}
