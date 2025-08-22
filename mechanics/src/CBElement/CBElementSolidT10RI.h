/*
 *  CBElementSolidT10RI.h
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 21.06.11.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#ifndef CB_ELEMENT_SOLID_T10_RI_H
#define CB_ELEMENT_SOLID_T10_RI_H

#include "CBElementSolidT10.h"

class CBElementSolidT10RI : public CBElementSolidT10
{
public:
    CBElementSolidT10RI(){}
    virtual ~CBElementSolidT10RI(){}
    static CBElement* New(){return(new CBElementSolidT10RI);}
    std::string GetType(){return(std::string("T10RI"));}
    virtual bool IsElementInverted();

protected:
private:
    CBStatus CalcNodalForcesHelperFunction(const TFloat* nodesCoords, const bool* boundaryConditions, const Matrix3<TFloat>* activeStress, TFloat* forces);

    typedef CBElement           Base;
    typedef CBElementSolidT10   Ancestor;

    CBElementSolidT10RI(const CBElementSolidT10RI &);
    void operator = (const CBElementSolidT10RI &);
};
#endif
