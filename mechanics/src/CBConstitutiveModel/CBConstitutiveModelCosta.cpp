/*
 *  CBConstitutiveModelCosta.cpp
 *  CardioMechanics
 *
 *  Created by Sonia Sassi on 07.04.2015.
 *  Copyright 2015 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#include "CBConstitutiveModelCosta.h"
#include "CBConstitutiveModelGuccione.h"

void CBConstitutiveModelCosta::Init(ParameterMap* parameters, TInt materialIndex)
{
    Base::Init(parameters,materialIndex);
    std::stringstream mi;
    mi << materialIndex;
    if(parameters->IsAvailable("Materials.Mat_" + mi.str() + ".Costa.C") == false && parameters->IsAvailable("Materials.Mat_Default.Costa.C") == true)
        c_ =  parameters->Get<double>("Materials.Mat_Default.Costa.C");
    else
        c_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".Costa.C");
    
    if(parameters->IsAvailable("Materials.Mat_" + mi.str() + ".Costa.bff") == false && parameters->IsAvailable("Materials.Mat_Default.Costa.bff") == true)
        bff_ =  parameters->Get<double>("Materials.Mat_Default.Costa.bff");
    else
        bff_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".Costa.bff");
    
    if(parameters->IsAvailable("Materials.Mat_" + mi.str() + ".Costa.bss") == false && parameters->IsAvailable("Materials.Mat_Default.Costa.bss") == true)
        bss_ =  parameters->Get<double>("Materials.Mat_Default.Costa.bss");
    else
        bss_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".Costa.bss");
    
    if(parameters->IsAvailable("Materials.Mat_" + mi.str() + ".Costa.bnn") == false && parameters->IsAvailable("Materials.Mat_Default.Costa.bnn") == true)
        bnn_ =  parameters->Get<double>("Materials.Mat_Default.Costa.bnn");
    else
        bnn_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".Costa.bnn");
    
    if(parameters->IsAvailable("Materials.Mat_" + mi.str() + ".Costa.bfs") == false && parameters->IsAvailable("Materials.Mat_Default.Costa.bfs") == true)
        bfs_ =  parameters->Get<double>("Materials.Mat_Default.Costa.bfs");
    else
        bfs_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".Costa.bfs");
    
    if(parameters->IsAvailable("Materials.Mat_" + mi.str() + ".Costa.bfn") == false && parameters->IsAvailable("Materials.Mat_Default.Costa.c5") == true)
        bfn_ =  parameters->Get<double>("Materials.Mat_Default.Costa.bfn");
    else
        bfn_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".Costa.bfn");
    
    if(parameters->IsAvailable("Materials.Mat_" + mi.str() + ".Costa.bns") == false && parameters->IsAvailable("Materials.Mat_Default.Costa.bns") == true)
        bns_ =  parameters->Get<double>("Materials.Mat_Default.Costa.cbns");
    else
        bns_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".Costa.bns");
    
    
    if(parameters->IsAvailable("Materials.Mat_" + mi.str() + ".Costa.K") == false && parameters->IsAvailable("Materials.Mat_Default.Costa.K") == true)
        k_ =  parameters->Get<double>("Materials.Mat_Default.Costa.K");
    else
        k_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".Costa.K");
}


CBStatus CBConstitutiveModelCosta::CalcEnergy(const Matrix3<TFloat>& deformationTensor, TFloat& energy)
{
    Matrix3<TFloat> greenStrainTensor = 0.5 * (deformationTensor * deformationTensor.GetTranspose() - identity_);
    TFloat*         e                 = greenStrainTensor.GetArray();
    TFloat          j                 = greenStrainTensor.Det();
    energy = c_ / 2.0 * (relaxedExp(bff_ * (e[0] * e[0]) + bss_*(e[4] * e[4]) + bnn_*(e[8] * e[8]) + 2*bfn_*(0.5*0.5*(e[2] + e[6])*(e[2] + e[6])) + 2*bfs_ *(0.5*0.5*(e[1]+ e[3])*(e[1]+ e[3])) + 2* bns_ *(0.5* 0.5* (e[7]+ e[5])* (e[7]+ e[5]))) - 1.0) + k_ / 2 * (j - 1) * (j - 1);
    if(std::isinf(energy))
        return(CBStatus::INFINITIVE);
    else if(deformationTensor.Det() <= 0)
        return(CBStatus::CORRUPT_ELEMENT);
    else
        return(CBStatus::SUCCESS);
}


CBStatus CBConstitutiveModelCosta::CalcPK2Stress(const Matrix3<TFloat>& deformationTensor, Matrix3<TFloat>& pk2Stress)
{
    
    if (!Base::ignoreCorruptElements_)
        if(deformationTensor.Det() <= 0)
            return(CBStatus::CORRUPT_ELEMENT);
    
    Matrix3<TFloat> rightGreenStrain  =  deformationTensor.GetTranspose() * deformationTensor;
    Matrix3<TFloat> greenStrainTensor = 0.5 * (rightGreenStrain - identity_);
    TFloat*         e                 = greenStrainTensor.GetArray();
    TFloat*         p                 = pk2Stress.GetArray();
    
    TFloat          expTerm = relaxedExp(bff_ * (e[0] * e[0]) + bss_*(e[4] * e[4]) + bnn_*(e[8] * e[8]) + 2*bfn_*(0.5*0.5*(e[2] + e[6])*(e[2] + e[6])) + 2*bfs_ *(0.5*0.5*(e[1]+ e[3])*(e[1]+ e[3])) + 2* bns_ *(0.5* 0.5* (e[7]+ e[5])* (e[7]+ e[5])));
    
    p[0] = c_ * bff_ * e[0] * expTerm;
    p[4] = c_ * bss_ * e[4] * expTerm;
    p[8] = c_ * bnn_ * e[8] * expTerm;
    
    p[1] = c_ * 0.5* bfs_ * (e[1] + e[3]) * expTerm;
    p[3] = c_ * 0.5* bfs_ * (e[1] + e[3]) * expTerm;
    
    p[2] = c_ * 0.5* bfn_ * (e[2] + e[6]) * expTerm;
    p[6] = c_ * 0.5* bfn_ * (e[2] + e[6]) * expTerm;
    
    p[5] = c_ * 0.5* bns_ * (e[5] + e[7]) * expTerm;
    p[7] = c_ * 0.5* bns_ * (e[5] + e[7]) * expTerm;
    
    
    
//#warning must be fixed !!! Highest Priority [ep901] WHAT MUST BE FIXED HERE???
    // Volume preservation
    
    Matrix3<TFloat> rightGreenStrainInverse = rightGreenStrain.GetInverse();
    TFloat          invariant3              = rightGreenStrain.Invariant3();
    //   pk2Stress += 2 * k_ * (log(invariant3)) * rightGreenStrainInverse; */
    
   pk2Stress += k_*(invariant3*invariant3 - 1) * rightGreenStrainInverse;
    
    if(std::isinf(expTerm))
        return(CBStatus::INFINITIVE);
    else if(std::isnan(expTerm))
        return(CBStatus::NOT_A_NUMBER);
    else
        return(CBStatus::SUCCESS);
}







