/*
 *  CBElementSolidT10RIT4.h
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 21.06.11.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

// Mixed element Taylor Hood tetrahedral element

#ifndef CB_ELEMENT_SOLID_T10_RI_T4_H
#define CB_ELEMENT_SOLID_T10_RI_T4_H

#include "CBElementSolidT10.h"

class CBElementSolidT10RIT4 : public CBElementSolidT10
{
public:
    CBElementSolidT10RIT4(){}
    virtual ~CBElementSolidT10RIT4(){}
    static CBElement* New();
    std::string GetType(){return(std::string("T10RIT4"));}
    unsigned int GetNumberOfNodesIndices(){if(isT10_)return(10);else return(4);}
    void UpdateShapeFunctions();
    bool IsElementInverted();
    
    void SwitchToT4(){isT10_ = false;}
    void SwitchToT10(){isT10_ = true;}
    
protected:
    TFloat dNdXt4_[12];
private:
    virtual CBStatus CalcNodalForcesHelperFunction(const TFloat* nodesCoords, const bool* boundaryConditions, const Matrix3<TFloat>* activeStress, TFloat* forces);
    void CalcT4ShapeFunctionsDerivatives();
    void CalcDeformationTensorsAtCentroidWithLocalBasis(const TFloat* nodesCoords, Matrix3<TFloat>& deformationTensors);
    
    typedef CBElement           Base;
    typedef CBElementSolidT10   Ancestor;
    
    CBElementSolidT10RIT4(const CBElementSolidT10RIT4 &);
    void operator=(const CBElementSolidT10RIT4 &);
    bool localStiffening_ = false;
    bool isProblemChild_ = false;
    bool isT10_ = true;
    
};

#endif
