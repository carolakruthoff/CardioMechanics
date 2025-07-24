/*! \file CBTensionEstimator.cpp
 * \brief implementation of CBTensionEstimator.cpp
 *
 * \ref CardioMechanics, CBTensionModel
 *
 * \author Ekaterina Kovacheva \n
 *  Copyright 2017 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 * \date  2019 18 September
 */


// includes
#include "CBTensionEstimator.h"
#include "DCCtrl.h"





//---------------------------------------------------------------------------------------------------
/*!
 Constructor
 \param e the finite element for which the tension model is initialized
 */

CBTensionEstimator::CBTensionEstimator( CBElement* e, ParameterMap* parameters )
{
    e_ = e;
    assert(e_ != nullptr);

    Tmax_ = e->GetMaterial()->GetProperties()->tensionMax_;
    activeTension_ = 0.0;

    // every element hat its own activeTensionForQuadPoints_ only for its QPs
    // there is no global activeTensionForQuadPoints_ vector
/*    activeTensionForQuadPoints_.resize( e_->GetNumberOfQuadraturePoints() );
    Tmax_ = e->GetMaterial()->GetProperties()->tensionMax_;
    int mi = e->GetMaterialIndex();
    coupledQP_ = parameters->Get<TInt>("Materials.Mat_" + std::to_string(mi) + ".ExternalTension.NumQPCoupling", 4);*/

}

//---------------------------------------------------------------------------------------------------
/*!
 Get the active stress vector of the quadrature points; length = nrOfElemetns * nrQPs
 \return active stress vector
 */
// std::vector<TFloat> CBTensionEstimator::GetActiveTension()
TFloat CBTensionEstimator::GetActiveTension()
{
    //ek717: added the multiplication with Tmax to be consistant with CalcActiveTension
    return activeTension_*Tmax_;
}


//---------------------------------------------------------------------------------------------------
/*!
 Calc Active Tension for the element based on the values asigned for the QPs ????? TODO
 return tension the value of the active tension for the current element
 */

TFloat CBTensionEstimator::CalcActiveTension(const math_pack::Matrix3<TFloat>& deformation, const TFloat time )
{ // Not actually used during calculations, but needed for the exporter – thus all QPs are avged
    // TFloat tension = 0.0;
    /*if (e_->GetType() == "T4") {
      tension = activeTensionForQuadPoints_[0];
    } else if (e_->GetType() == "T10"){
      for (int i = 0; i < coupledQP_; i++){
      tension += activeTensionForQuadPoints_[i]/coupledQP_;
      }
    }
*/
    return activeTension_*Tmax_;
}


CBStatus CBTensionEstimator::SetActiveTensionAtQuadraturePoint( TInt indexQP, TFloat activeTension ){
    
    CBStatus rc = CBStatus::SUCCESS;
    
    activeTension_ = activeTension;
    
    return rc;
}



