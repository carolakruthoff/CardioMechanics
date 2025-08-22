//
//  CBElementSolidT10T4.cpp
//  CardioMechanics_unstable
//
//  Created by Thomas Fritz on 16.01.14.
//
//

#include "CBElementSolidT10T4.h"
#include "CBElementAdapter.h"
#include "DCCtrl.h"

CBElement* CBElementSolidT10T4::New()
{
    return(new CBElementSolidT10T4);
}


bool CBElementSolidT10T4::IsElementInverted()
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

CBStatus CBElementSolidT10T4::CalcNodalForcesHelperFunction(const TFloat* nodesCoords, const bool* boundaryConditions, const Matrix3<TFloat>* activeStress, TFloat* forces)
{
   
    CBStatus rc;
    Matrix3<TFloat> deformationTensors[5];
    Matrix3<TFloat> pk2Stress[5];
    
    CalcDeformationTensorsAtQuadraturePointsWithLocalBasis(nodesCoords, deformationTensors);
    CalcDeformationTensorsAtCentroidWithLocalBasisWithT4ShapeFunctions(nodesCoords, deformationTensors[4]);
    
    for(int i = 0; i < 5; i++)
    {
        rc = Base::material_->GetConstitutiveModel()->CalcPK2Stress(deformationTensors[i], pk2Stress[i]);
        
//        if(!isDefect_)
//        {
//            if(rc == CBStatus::CORRUPT_ELEMENT)
//            {
//                isDefect_ = true;
//                std::cout << "Corrupt element " << GetIndex() << "\n";
//                return CBStatus::SUCCESS;
//            }
//            if(rc == CBStatus::CRITICAL_STRESS)
//            {
//                isDefect_ = true;
//                std::cout << "Critical Stress in element " << GetIndex() << "\n";
//                return CBStatus::SUCCESS;
//            }
//        }
//        else
//            return CBStatus::SUCCESS;
        
       
        pk2Stress[i] = pk2Stress[i] + *activeStress;
        pk2Stress[i] = pk2Stress[i] * deformationTensors[i].GetTranspose();
        
        if(i!=4)
            pk2Stress[i] = GetBasisAtQuadraturePoint(i+1)->GetInverse().GetTranspose() * pk2Stress[i] * GetBasisAtQuadraturePoint(i+1)->GetTranspose();
        else
            pk2Stress[4] = GetBasisAtQuadraturePoint(0)->GetInverse().GetTranspose() * pk2Stress[i] * GetBasisAtQuadraturePoint(0)->GetTranspose();
        
        
    }
  
    
    TFloat  t[3];
    // Calculate forces;
    
    for(int i=0; i <30; i++)
        forces[i] = 0;
    
    if(isT10_)
        for(int i = 0; i < 10; i++)
        {
            t[0] = 0;
            t[1] = 0;
            t[2] = 0;
            
            for(int j = 0; j < 4; j++)
            {
                t[0] += 0.25 * (detJ_/6.0) * (dNdXW_[30*j+3*i+0] * pk2Stress[j].Get(0, 0)  + dNdXW_[30*j+3*i+1] * pk2Stress[j].Get(1, 0) + dNdXW_[30*j+3*i+2] * pk2Stress[j].Get(2, 0));
                t[1] += 0.25 * (detJ_/6.0) * (dNdXW_[30*j+3*i+0] * pk2Stress[j].Get(0, 1)  + dNdXW_[30*j+3*i+1] * pk2Stress[j].Get(1, 1) + dNdXW_[30*j+3*i+2] * pk2Stress[j].Get(2, 1));
                t[2] += 0.25 * (detJ_/6.0) * (dNdXW_[30*j+3*i+0] * pk2Stress[j].Get(0, 2)  + dNdXW_[30*j+3*i+1] * pk2Stress[j].Get(1, 2) + dNdXW_[30*j+3*i+2] * pk2Stress[j].Get(2, 2));
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
    else
        for(int i = 0; i < 4; i++)
        {
            if(boundaryConditions[3*i] == 0)
                forces[3*i] = Ancestor::initialVolume_*(pk2Stress[4].Get(0, 0) *  dNdXt4_[3*i+0] + pk2Stress[4].Get(1, 0) *  dNdXt4_[3*i+1] + pk2Stress[4].Get(2, 0) *  dNdXt4_[3*i+2]);
            else
                forces[3*i] = 0;
            
            if(boundaryConditions[3*i+1] == 0)
                forces[3*i+1] = Ancestor::initialVolume_*(pk2Stress[4].Get(0, 1) *  dNdXt4_[3*i+0] + pk2Stress[4].Get(1, 1) *  dNdXt4_[3*i+1] + pk2Stress[4].Get(2, 1) *  dNdXt4_[3*i+2]);
            else
                forces[3*i+1] = 0;
            
            if(boundaryConditions[3*i+2] == 0)
                forces[3*i+2] = Ancestor::initialVolume_*(pk2Stress[4].Get(0, 2) *  dNdXt4_[3*i+0] + pk2Stress[4].Get(1, 2) *  dNdXt4_[3*i+1] + pk2Stress[4].Get(2, 2) *  dNdXt4_[3*i+2]);
            else
                forces[3*i+2] = 0;
        }

    return CBStatus::SUCCESS;
}