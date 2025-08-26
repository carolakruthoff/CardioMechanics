/* -------------------------------------------------------

   CBTensionModelExternal.cpp

   Ver. 1.0.0

   Created:       Ekaterina Kovacheva (25.04.2018)
   Last modified: Tobias Gerach (21.12.2022)

   Institute of Biomedical Engineering
   Karlsruhe Institute of Technology (KIT)

   http://www.ibt.kit.edu

   Copyright 2000-2009 - All rights reserved.

   ------------------------------------------------------ */


// includes
#include "CBTensionModelExternal.h"
#include "DCCtrl.h"


// ---------------------------------------------------------------------------------------------------
/*!
   Constructor
   \param e the finite element for which the tension model is initialized
 */

CBTensionModelExternal::CBTensionModelExternal(CBElement *e, ParameterMap *parameters) {
  e_ = e;
  assert(e_ != nullptr);

  // every element has its own activeTensionForQuadPoints_ only for its QPs
  // there is no global activeTensionForQuadPoints_ vector
  activeTensionForQuadPoints_.resize(e_->GetNumberOfQuadraturePoints() );
  Tmax_ = e->GetMaterial()->GetProperties()->tensionMax_;
}

// ---------------------------------------------------------------------------------------------------

/*!
   Set Active Tension for one Quadrature Point (QP)
   \param indexQP the index of the QP
   \param activeTension the value of the active tension
   \return status
 */
CBStatus CBTensionModelExternal::SetActiveTensionAtQuadraturePoint(TInt indexQP, TFloat activeTension) {
  CBStatus rc = CBStatus::SUCCESS;

  if (e_->GetType() == "T4") {
    if (indexQP == 0) {
      activeTensionForQuadPoints_[indexQP] = activeTension;
    } else {
      throw std::runtime_error(
              "CBTensionAcCELLerateModel::SetActiveTensionForOneQuaraturePoint: try to access a not existing qudrature point for T4 elements.");
    }
  } else if (e_->GetType() == "T10") {
    if (indexQP < 5) {
      activeTensionForQuadPoints_[indexQP] = activeTension;
    } else {
      throw std::runtime_error(
              "CBTensionAcCELLerateModel::SetActiveTensionForOneQuaraturePoint: try to access a not existing qudrature point for T10 elements.");
    }
  } else {
    rc = CBStatus::FAILED;
  }

  return rc;
} // CBTensionModelExternal::SetActiveTensionAtQuadraturePoint

// ---------------------------------------------------------------------------------------------------

/*!
   Calc Active Tension for the element based on the values asigned for the QPs ????? TODO
   return tension the value of the active tension for the current element
 */
TFloat CBTensionModelExternal::CalcActiveTension(const math_pack::Matrix3<TFloat> &deformation,
                                                 const TFloat time) { // Not actually used during calculations, but needed for the exporter – thus all QPs are avged
  TFloat tension = 0.0;

  if (e_->GetType() == "T4") {
    tension = activeTensionForQuadPoints_[0];
  } else if (e_->GetType() == "T10") {
    for (int i = 1; i < 5; i++) {
      tension += activeTensionForQuadPoints_[i]/4;
    }
  }

  return tension*Tmax_;
}

// ---------------------------------------------------------------------------------------------------

/*!
   Get the active stress vector of the quadrature points; length = nrOfElemetns * nrQPs
   \return active stress vector
 */
std::vector<TFloat> CBTensionModelExternal::GetActiveTensionVectorOfQPs() {
  return activeTensionForQuadPoints_;
}

// ---------------------------------------------------------------------------------------------------

/*!
   //!< Active stress for one quadratire point
   //return activeStress 3x3 matrix = | activeTension  0  0 |
 |       0        0  0 |
 |       0        0  0 |

 */
math_pack::Matrix3<TFloat>  CBTensionModelExternal::CalcActiveStressAtQuadraturePoint(
  const math_pack::Matrix3<TFloat> &deformation, const TFloat time, TInt indexQP) {
  TFloat activeTensionQP = activeTensionForQuadPoints_[indexQP]*Tmax_;

  return math_pack::Matrix3<TFloat>(activeTensionQP * fibreRatio_.Get(0), 0, 0,
                                    0, activeTensionQP * fibreRatio_.Get(1), 0,
                                    0, 0, activeTensionQP * fibreRatio_.Get(2));
}
