/* -------------------------------------------------------

   CBTensionModelExternal.h

   Ver. 1.0.0

   Created:       Ekaterina Kovacheva (25.04.2018)
   Last modified: Tobias Gerach (21.12.2022)

   Institute of Biomedical Engineering
   Karlsruhe Institute of Technology (KIT)

   http://www.ibt.kit.edu

   Copyright 2000-2009 - All rights reserved.

   ------------------------------------------------------ */

#pragma once

// includes
#include "CBTensionModel.h"

#include "CBDataCtrl.h"


class CBTensionModelExternal : public CBTensionModel {
 protected:
  CBElement *e_;  //!< solid element object

  std::vector<TFloat> activeTensionForQuadPoints_;

  Vector3<TFloat> fibreRatio_ = {1.0, 0.0, 0.0};

 public:
  //!< Constructor CBTensionAcCELLerateModel
  CBTensionModelExternal(CBElement *, ParameterMap *);

  //!< Destructor CBTensionAcCELLerateModel
  CBTensionModelExternal() {}

  //!< Get the active stress vector of the quadrature points
  std::vector<TFloat> GetActiveTensionVectorOfQPs();

  //!< Set Active Tension for one Quadrature Point
  CBStatus SetActiveTensionAtQuadraturePoint(TInt, TFloat) override;

  //!< Active tension computation
  virtual TFloat CalcActiveTension(const math_pack::Matrix3<TFloat> &, const TFloat) override;

  //!< Active stress for one quadrature point
  virtual math_pack::Matrix3<TFloat>  CalcActiveStressAtQuadraturePoint(const math_pack::Matrix3<TFloat>&, const TFloat,
                                                                        TInt) override;

  //! Set Fibre Ratio
  virtual void SetfibreRatio(Vector3<TFloat> ffRatio) override {fibreRatio_ = ffRatio;}
};
