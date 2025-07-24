/*! \file CBTensionEstimator.h
 * \brief CBTensionEstimator to estimate tension with the inverse problem
 *
 * \ref CardioMechanics, CBTensionModel
 *
 * \author Ekaterina Kovacheva \n
 *  Copyright 2017 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 * \date  2019 18 September
 */

#pragma once

// includes
#include "CBTensionModel.h"

#include "CBDataCtrl.h"


class CBTensionEstimator : public CBTensionModel {
    
protected:
    
    CBElement* e_;  //!< solid element object
    TFloat activeTension_;

public:

    //!< Constructor CBTensionAcCELLerateModel
    CBTensionEstimator( CBElement*, ParameterMap*);
    
    //!< Destructor CBTensionAcCELLerateModel
    CBTensionEstimator() { }

    //!< Get the active stress vector of the quadrature points
    // std::vector<TFloat> GetActiveTension();
    TFloat GetActiveTension();

    //!< Active tension computation
    virtual TFloat CalcActiveTension(const math_pack::Matrix3<TFloat>& , const TFloat) override;
    
    CBStatus SetActiveTensionAtQuadraturePoint( int indexQP, TFloat activeTension) override;
    

};
