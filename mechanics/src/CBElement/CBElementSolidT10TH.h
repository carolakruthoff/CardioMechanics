/*
 *  CBElementSolidT10TH.h
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 21.06.11.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

// Mixed element Taylor Hood tetrahedral element

#ifndef CB_ELEMENT_SOLID_T10_TH_H
#define CB_ELEMENT_SOLID_T10_TH_H

#include "CBElementSolidT10.h"

class CBElementSolidT10TH : public CBElementSolidT10
{
public:
    CBElementSolidT10TH(){}
    virtual ~CBElementSolidT10TH(){}
    static CBElement* New(){return(new CBElementSolidT10TH);}
    std::string GetType(){return(std::string("T10TH"));}
    
protected:
private:
    virtual CBStatus CalcNodalForcesHelperFunction(const TFloat* nodesCoords, const bool* boundaryConditions, const Matrix3<TFloat>* activeStress, TFloat* forces);

    typedef CBElement           Base;
    typedef CBElementSolidT10   Ancestor;

    CBElementSolidT10TH(const CBElementSolidT10TH &);
    void operator=(const CBElementSolidT10TH &);
};

#endif
