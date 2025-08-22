/*
 *  CBElementSolidT10T4.h
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 21.06.11.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

// Mixed element Taylor Hood tetrahedral element

#ifndef CB_ELEMENT_SOLID_T4_T10_H
#define CB_ELEMENT_SOLID_T4_T10_H

#include "CBElementSolidT10.h"

class CBElementSolidT10T4 : public CBElementSolidT10
{
public:
    CBElementSolidT10T4() : CBElementSolidT10() {}
    virtual ~CBElementSolidT10T4(){}
    static CBElement* New();
    std::string GetType(){return(std::string("T10T4"));}
    unsigned int GetNumberOfNodesIndices(){if(isT10_)return(10);else return(4);}
    bool IsElementInverted();
    
    void SwitchToT4(){isT10_ = false;}
    void SwitchToT10(){isT10_ = true;}
    
protected:
private:
    virtual CBStatus CalcNodalForcesHelperFunction(const TFloat* nodesCoords, const bool* boundaryConditions, const Matrix3<TFloat>* activeStress, TFloat* forces);
    
    typedef CBElement           Base;
    typedef CBElementSolidT10   Ancestor;
    
    CBElementSolidT10T4(const CBElementSolidT10T4 &);
    void operator=(const CBElementSolidT10T4 &);
    bool isT10_ = true;
    
};

#endif
