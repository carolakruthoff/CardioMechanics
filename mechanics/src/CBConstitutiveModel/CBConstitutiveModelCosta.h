/*
 *  CBConstitutiveModelCosta.h
 *  CardioMechanics
 *
 *  Created by Sonia Sassi on 07.04.2015.
 *  Copyright 2015 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

//    "Modelling cardiac mechanical properties in three dimensions" (2001)
//    Costa, K.D. and Holmes, J.W. and McCulloch, A.D.
//    Philosophical Transactions of the Royal Society of London. Series A: Mathematical, Physical and Engineering Sciences Vol 359,num 1783,page 1233,
//
//    "Mechanische Modellierung und Kalibrierung von Muskelfaser-Schichten im menschlichen Herzen" (2015), Sonia Sassi, Bachelor Thesis


    


#ifndef CB_CONSTITUTIVE_MODEL_COSTA
#define CB_CONSTITUTIVE_MODEL_COSTA

#include<string>
#include <sstream>

#include "Matrix3.h"

#include "ParameterMap.h"

#include "CBConstitutiveModel.h"



class CBConstitutiveModelCosta : public CBConstitutiveModel
{
public:
    void Init(ParameterMap* parameters, TInt materialIndex);
    CBConstitutiveModelCosta() : c_(0), bff_(0), bss_(0), bnn_(0), bfs_(0), bfn_(0), bns_(0), k_(0) {identity_.SetToIdentityMatrix(); }
    CBStatus CalcEnergy(const Matrix3<TFloat>& deformationTensor, TFloat& energy);
    CBStatus CalcPK2Stress(const Matrix3<TFloat>& deformationTensor, Matrix3<TFloat>& pk2Stress);
protected:
private:
    typedef CBConstitutiveModel Base;
    TFloat          c_, bff_,bss_,bnn_,bfs_,bfn_,bns_, k_;
    Matrix3<TFloat> identity_;
};

#endif
