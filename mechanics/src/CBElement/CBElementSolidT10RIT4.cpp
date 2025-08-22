//
//  CBElementSolidT10RIT4.cpp
//  CardioMechanics_unstable
//
//  Created by Thomas Fritz on 16.01.14.
//
//

#include "CBElementSolidT10RIT4.h"
#include "CBElementAdapter.h"

CBElement* CBElementSolidT10RIT4::New()
{
    return(new CBElementSolidT10RIT4);
}

void CBElementSolidT10RIT4::CalcT4ShapeFunctionsDerivatives()
{
    
    TFloat nodesCoords[30];
    TInt   nodesCoordsIndices[30];
    
    
    GetNodesCoordsIndices(nodesCoordsIndices);
    Base::adapter_->GetNodesCoords(30, nodesCoordsIndices, nodesCoords);
    
    TFloat z43 = (nodesCoords[11] - nodesCoords[8]);
    TFloat z42 = (nodesCoords[11] - nodesCoords[5]);
    TFloat z41 = (nodesCoords[11] - nodesCoords[2]);
    TFloat z32 = (nodesCoords[ 8] - nodesCoords[5]);
    TFloat z31 = (nodesCoords[ 8] - nodesCoords[2]);
    TFloat z21 = (nodesCoords[ 5] - nodesCoords[2]);
    
    TFloat y43 = (nodesCoords[10] - nodesCoords[7]);
    TFloat y42 = (nodesCoords[10] - nodesCoords[4]);
    TFloat y41 = (nodesCoords[10] - nodesCoords[1]);
    TFloat y32 = (nodesCoords[ 7] - nodesCoords[4]);
    TFloat y31 = (nodesCoords[ 7] - nodesCoords[1]);
    TFloat y21 = (nodesCoords[ 4] - nodesCoords[1]);
    
    
    TFloat x41 = (nodesCoords[ 9] - nodesCoords[0]);
    TFloat x31 = (nodesCoords[ 6] - nodesCoords[0]);
    TFloat x21 = (nodesCoords[ 3] - nodesCoords[0]);
    
    TFloat v = x21 * (y31 * z41 - y41 * z31) + y21 * (x41 * z31 - x31 * z41) + z21 * (x31 * y41 - x41 * y31);
    
    dNdXt4_[ 0] = 1.0 /v*( nodesCoords[4] * z43 - nodesCoords[7] * z42 + nodesCoords[10] * z32);
    dNdXt4_[ 3] = 1.0 /v*(-nodesCoords[1] * z43 + nodesCoords[7] * z41 - nodesCoords[10] * z31);
    dNdXt4_[ 6] = 1.0 /v*( nodesCoords[1] * z42 - nodesCoords[4] * z41 + nodesCoords[10] * z21);
    dNdXt4_[ 9] = 1.0 /v*(-nodesCoords[1] * z32 + nodesCoords[4] * z31 - nodesCoords[ 7] * z21);
    
    dNdXt4_[ 1] = 1.0 /v*(-nodesCoords[3] * z43 + nodesCoords[6] * z42 - nodesCoords[ 9] * z32);
    dNdXt4_[ 4] = 1.0 /v*( nodesCoords[0] * z43 - nodesCoords[6] * z41 + nodesCoords[ 9] * z31);
    dNdXt4_[ 7] = 1.0 /v*(-nodesCoords[0] * z42 + nodesCoords[3] * z41 - nodesCoords[ 9] * z21);
    dNdXt4_[10] = 1.0 /v*( nodesCoords[0] * z32 - nodesCoords[3] * z31 + nodesCoords[ 6] * z21);
    
    dNdXt4_[ 2] = 1.0 /v*( nodesCoords[3] * y43 - nodesCoords[6] * y42 + nodesCoords[ 9] * y32);
    dNdXt4_[ 5] = 1.0 /v*(-nodesCoords[0] * y43 + nodesCoords[6] * y41 - nodesCoords[ 9] * y31);
    dNdXt4_[ 8] = 1.0 /v*( nodesCoords[0] * y42 - nodesCoords[3] * y41 + nodesCoords[ 9] * y21);
    dNdXt4_[11] = 1.0 /v*(-nodesCoords[0] * y32 + nodesCoords[3] * y31 - nodesCoords[ 6] * y21);
}


void CBElementSolidT10RIT4::UpdateShapeFunctions()
{
    CBElementSolidT10::UpdateShapeFunctions();
    CalcT4ShapeFunctionsDerivatives();
}

void CBElementSolidT10RIT4::CalcDeformationTensorsAtCentroidWithLocalBasis(const TFloat* nodesCoords, Matrix3<TFloat>& deformationTensor)
{
    
    TFloat* f = deformationTensor.GetArray();
    
    for(unsigned int i = 0; i < 3; i++)
    {
        f[i]   = dNdXt4_[i] * nodesCoords[0] + dNdXt4_[i+3] * nodesCoords[3] + dNdXt4_[i+6] * nodesCoords[6] + dNdXt4_[i+9] * nodesCoords[ 9];
        f[3+i] = dNdXt4_[i] * nodesCoords[1] + dNdXt4_[i+3] * nodesCoords[4] + dNdXt4_[i+6] * nodesCoords[7] + dNdXt4_[i+9] * nodesCoords[10];
        f[6+i] = dNdXt4_[i] * nodesCoords[2] + dNdXt4_[i+3] * nodesCoords[5] + dNdXt4_[i+6] * nodesCoords[8] + dNdXt4_[i+9] * nodesCoords[11];
    }
    
    deformationTensor = GetBasisAtQuadraturePoint(0)->GetTranspose() * deformationTensor * GetBasisAtQuadraturePoint(0)->GetInverse().GetTranspose(); //
    
}

bool CBElementSolidT10RIT4::IsElementInverted()
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

CBStatus CBElementSolidT10RIT4::CalcNodalForcesHelperFunction(const TFloat* nodesCoords, const bool* boundaryConditions, const Matrix3<TFloat>* activeStress, TFloat* forces)
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
    if(isT10_)
    {
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
    }
    else
    {
        Matrix3<TFloat> pk2Stress = pk2Deviatoric[4]+pk2Hydrostatic[4];
        for(int i = 0; i < 4; i++)
        {
            if(boundaryConditions[3*i] == 0)
                forces[3*i] = Ancestor::initialVolume_*(pk2Stress.Get(0, 0) *  dNdXt4_[3*i+0] + pk2Stress.Get(1, 0) *  dNdXt4_[3*i+1] + pk2Stress.Get(2, 0) *  dNdXt4_[3*i+2]);
            else
                forces[3*i] = 0;
            
            if(boundaryConditions[3*i+1] == 0)
                forces[3*i+1] = Ancestor::initialVolume_*(pk2Stress.Get(0, 1) *  dNdXt4_[3*i+0] + pk2Stress.Get(1, 1) *  dNdXt4_[3*i+1] + pk2Stress.Get(2, 1) *  dNdXt4_[3*i+2]);
            else
                forces[3*i+1] = 0;
            
            if(boundaryConditions[3*i+2] == 0)
                forces[3*i+2] = Ancestor::initialVolume_*(pk2Stress.Get(0, 2) *  dNdXt4_[3*i+0] + pk2Stress.Get(1, 2) *  dNdXt4_[3*i+1] + pk2Stress.Get(2, 2) *  dNdXt4_[3*i+2]);
            else
                forces[3*i+2] = 0;
        }
    }
    
    
    return(rc);
}

