/*
 *  CBElementSolidT10RI.cpp
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 13.09.11.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#include "CBElementSolidT10RI.h"
#include "CBElementAdapter.h"

bool CBElementSolidT10RI::IsElementInverted()
{
    
    TFloat          nodesCoords[30];
    TInt            nodesCoordsIndices[30];
    Matrix3<TFloat> deformationTensors[5];
    
    GetNodesCoordsIndices(nodesCoordsIndices);
    Base::adapter_->GetNodesCoords(30, nodesCoordsIndices, nodesCoords);
    CalcDeformationTensorsAtQuadraturePointsWithLocalBasis(nodesCoords, deformationTensors);
    CalcDeformationTensorsAtCentroidWithLocalBasis(nodesCoords, deformationTensors[4]);
    
    for(int i=0; i<5; i++)
        if(deformationTensors[i].Det() < 0)
            return true;
    
    return false;
}

CBStatus CBElementSolidT10RI::CalcNodalForcesHelperFunction(const TFloat* nodesCoords, const bool* boundaryConditions, const Matrix3<TFloat>* activeStress,  TFloat* forces)
{
    CBStatus rc;

    Matrix3<TFloat> deformationTensors[5];
    Matrix3<TFloat> pk2Deviatoric[5];
    Matrix3<TFloat> pk2Hydrostatic[5];

    CalcDeformationTensorsAtQuadraturePointsWithLocalBasis(nodesCoords, deformationTensors);
    CalcDeformationTensorsAtCentroidWithLocalBasis(nodesCoords, deformationTensors[4]);
    
    for(int i = 0; i < 5; i++)
    {
        Matrix3<TFloat> pk2Stress;
        rc = Base::material_->GetConstitutiveModel()->CalcPK2Stress(deformationTensors[i], pk2Stress);
        if(rc != CBStatus::SUCCESS)
            return(rc);
        
        pk2Stress += *activeStress;
        
        pk2Hydrostatic[i] = 1.0/3.0 * pk2Stress.Trace() * Matrix3<TFloat>(1,0,0,0,1,0,0,0,1);
        pk2Deviatoric[i] = pk2Stress - pk2Hydrostatic[i];
        
 //       CBElementSolid::RepairDeformationTensorIfInverted(deformationTensors[i]);
        
        pk2Deviatoric[i] = pk2Deviatoric[i] * deformationTensors[i].GetTranspose();
        pk2Hydrostatic[i] = pk2Hydrostatic[i] * deformationTensors[i].GetTranspose();
        
        if(i!=4)
            pk2Deviatoric[i] = GetBasisAtQuadraturePoint(i+1)->GetInverse().GetTranspose() * pk2Deviatoric[i] * GetBasisAtQuadraturePoint(i+1)->GetTranspose();
        else
            pk2Hydrostatic[4] = GetBasisAtQuadraturePoint(0)->GetInverse().GetTranspose() * pk2Hydrostatic[i] * GetBasisAtQuadraturePoint(0)->GetTranspose();
        
    }
    
    TFloat  t[3];
    // Calculate forces;
    for(int i = 0; i < 10; i++)
    {
        t[0] = 0.25 * (Ancestor::detJ_/6.0) * (Ancestor::dNdXW_[120+3*i+0] *  pk2Hydrostatic[4].Get(0, 0) + Ancestor::dNdXW_[120+3*i+1] * pk2Hydrostatic[4].Get(1, 0) + Ancestor::dNdXW_[120+3*i+2] * pk2Hydrostatic[4].Get(2, 0));
        t[1] = 0.25 * (Ancestor::detJ_/6.0) * (Ancestor::dNdXW_[120+3*i+0] *  pk2Hydrostatic[4].Get(0, 1) + Ancestor::dNdXW_[120+3*i+1] * pk2Hydrostatic[4].Get(1, 1) + Ancestor::dNdXW_[120+3*i+2] * pk2Hydrostatic[4].Get(2, 1));
        t[2] = 0.25 * (Ancestor::detJ_/6.0) * (Ancestor::dNdXW_[120+3*i+0] *  pk2Hydrostatic[4].Get(0, 2) + Ancestor::dNdXW_[120+3*i+1] * pk2Hydrostatic[4].Get(1, 2) + Ancestor::dNdXW_[120+3*i+2] * pk2Hydrostatic[4].Get(2, 2));
       
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

    return(rc);
}
