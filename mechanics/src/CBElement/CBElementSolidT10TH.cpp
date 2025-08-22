/*
 *  CBElementSolidT10TH.cpp
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 13.09.11.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#include "CBElementSolidT10TH.h"
#include "CBElementAdapter.h"

CBStatus CBElementSolidT10TH::CalcNodalForcesHelperFunction(const TFloat* nodesCoords, const bool* boundaryConditions, const Matrix3<TFloat>* activeStress, TFloat* forces)
{
    
    CBStatus rc;
    Matrix3<TFloat> deformationTensors[5];
    
    Matrix3<TFloat> d1[5];
    Matrix3<TFloat> d2[5];
    
    Matrix3<TFloat> pk2Deviatoric[5];
    Matrix3<TFloat> pk2Hydrostatic[5];
    
    CalcDeformationTensorsAtQuadraturePointsWithLocalBasis(nodesCoords, deformationTensors);
    CalcDeformationTensorsAtCentroidWithLocalBasisWithT4ShapeFunctions(nodesCoords, deformationTensors[4]);
    
    
    for(int i = 0; i < 5; i++)
    {
        Matrix3<TFloat> pk2Stress;
   //     CBElementSolid::RepairDeformationTensorIfInverted(deformationTensors[i]);

        d1[i] = 1.0/3.0 * deformationTensors[i].Trace() * Matrix3<TFloat>(1,0,0,0,1,0,0,0,1);
        d2[i] = deformationTensors[i] - d1[i];
        
        rc = Base::material_->GetConstitutiveModel()->CalcPK2Stress(d1[i], pk2Hydrostatic[i]);
        if(rc != CBStatus::SUCCESS)
            return(rc);
        
        rc = Base::material_->GetConstitutiveModel()->CalcPK2Stress(d2[i], pk2Deviatoric[i]);
        if(rc != CBStatus::SUCCESS)
            return(rc);
        
        pk2Deviatoric[i] = pk2Deviatoric[i] * d1[i].GetTranspose();
        pk2Hydrostatic[i] = pk2Hydrostatic[i] * d2[i].GetTranspose();
        
        pk2Deviatoric[i] += *activeStress;
        
        if(i!=4)
        {
            pk2Deviatoric[i] = GetBasisAtQuadraturePoint(i+1)->GetInverse().GetTranspose() * pk2Deviatoric[i] * GetBasisAtQuadraturePoint(i+1)->GetTranspose();
            pk2Hydrostatic[i] = GetBasisAtQuadraturePoint(i+1)->GetInverse().GetTranspose() * pk2Hydrostatic[i] * GetBasisAtQuadraturePoint(i+1)->GetTranspose();
        }
        else
        {
            pk2Deviatoric[4] = GetBasisAtQuadraturePoint(0)->GetInverse().GetTranspose() * pk2Deviatoric[4] * GetBasisAtQuadraturePoint(0)->GetTranspose();
            pk2Hydrostatic[4] = GetBasisAtQuadraturePoint(0)->GetInverse().GetTranspose() * pk2Hydrostatic[4] * GetBasisAtQuadraturePoint(0)->GetTranspose();
        }
        
    }
    
    TFloat  t[3];
    // Calculate forces;
    for(int i = 0; i < 10; i++)
    {
        if(i < 4)
        {
            t[0] = Ancestor::initialVolume_*(pk2Hydrostatic[4].Get(0,0) * dNdXt4_[3*i+0] + pk2Hydrostatic[4].Get(1,0) * dNdXt4_[3*i+1] + pk2Hydrostatic[4].Get(2,0) * dNdXt4_[3*i+2]);
            t[1] = Ancestor::initialVolume_*(pk2Hydrostatic[4].Get(0,1) * dNdXt4_[3*i+0] + pk2Hydrostatic[4].Get(1,1) * dNdXt4_[3*i+1] + pk2Hydrostatic[4].Get(2,1) * dNdXt4_[3*i+2]);
            t[2] = Ancestor::initialVolume_*(pk2Hydrostatic[4].Get(0,2) * dNdXt4_[3*i+0] + pk2Hydrostatic[4].Get(1,2) * dNdXt4_[3*i+1] + pk2Hydrostatic[4].Get(2,2) * dNdXt4_[3*i+2]);
        }
        else
        {
            t[0] = 0;
            t[1] = 0;
            t[2] = 0;
        }
        
        for(int j = 0; j < 4; j++)
        {
            t[0] += 0.25 * (Ancestor::detJ_/6.0) * (Ancestor::dNdXW_[30*j+3*i+0] *  pk2Deviatoric[j].Get(0, 0) + Ancestor::dNdXW_[30*j+3*i+1] * pk2Deviatoric[j].Get(1, 0) + Ancestor::dNdXW_[30*j+3*i+2] * pk2Deviatoric[j].Get(2, 0));
            t[1] += 0.25 * (Ancestor::detJ_/6.0) * (Ancestor::dNdXW_[30*j+3*i+0] *  pk2Deviatoric[j].Get(0, 1) + Ancestor::dNdXW_[30*j+3*i+1] * pk2Deviatoric[j].Get(1, 1) + Ancestor::dNdXW_[30*j+3*i+2] * pk2Deviatoric[j].Get(2, 1));
            t[2] += 0.25 * (Ancestor::detJ_/6.0) * (Ancestor::dNdXW_[30*j+3*i+0] *  pk2Deviatoric[j].Get(0, 2) + Ancestor::dNdXW_[30*j+3*i+1] * pk2Deviatoric[j].Get(1, 2) + Ancestor::dNdXW_[30*j+3*i+2] * pk2Deviatoric[j].Get(2, 2));
        }
        
        if(boundaryConditions[3*i] == 0)
            forces[3*i  ] = t[0];
        else
            forces[3*i  ] = 0;
        
        if(boundaryConditions[3*i+1] == 0)
            forces[3*i+1] = t[1];
        else
            forces[3*i+1] = 0;
        
        if(boundaryConditions[3*i+2] == 0)
            forces[3*i+2] = t[2];
        else
            forces[3*i+2] = 0;
    }

    return (rc);
}



