/*
 * ; *  CBWindkessel3.h
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 12.05.11.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#ifndef CB_WINDKESSEL_3
#define CB_WINDKESSEL_3

#include "Matrix3.h"

#include "CBElementCavity.h"
#include "CBCirculatorySystemModel.h"
#include "ParameterMap.h"

#include <string>

/*
 * See:
 *
 * Total arterial inertance as the fourth element of the
 * windkessel model
 * Nikos Stergiopulos, Berend E. Westerhof and Nico Westerhof
 * Am J Physiol Heart Circ Physiol 276:H81-H88, 1999.
 */

using namespace math_pack;

class CBWindkessel3 : public CBCirculatorySystemModel
{
public:
    CBWindkessel3(CBCavity* cavity, ParameterMap* parameters);
    TFloat Update(TFloat time, TFloat vol) override;

    TFloat GetPressure() override {return(ventricularPressure_);}
    TFloat CalcPressure(TFloat vol) override;
    TFloat CalcPressureDerivative(TFloat vol) override;
    void Prepare() override;

    TFloat GetPreparationProgress() override {return(time_/(preloadTime1_+preloadTime2_));}
    void StepBack() override;
    void WriteToFile(TFloat time) override;
    void InitInitialVolume() override;

protected:
private:
    TFloat            ventricularPressure_ = 0;
    TFloat            volumePressureWork_ = 0;
    TFloat            lastVolumePressureWork_ = 0;
    TFloat            isovolumetricRelaxationCavityVolume_ = 0;
    TFloat            atrialPressure_             = 0;
    TFloat            aorticPressure_             = 0;
    TFloat            beatInterval_               = 0;
    TFloat            volume_                     = 0;
    TFloat            initialVolume_              = 0;
    TFloat            initialAorticPressure_      = 0;
    TFloat            initialVentricularPressure_ = 0;
    TFloat            lastVentricularPressure_    = 0;
    TFloat            lastAorticPressure_         = 0;
    TFloat            lastVolume_                 = 0;
    TFloat            preloadTime1_               = 0;
    TFloat            preloadTime2_               = 0;
    TFloat            prevVolume_                 = 0;
    TFloat            lastPrevVolume_             = 0;
    TFloat            bulkModulus_                = 0;
    TFloat            dt_                         = 0;
    TFloat            prevDt_                     = 0;
    TFloat            lastPrevDt_                 = 0;
    TFloat            minSystolicPressure_        = 0;
    TFloat            maxVentricularPressure_     = 0;
    TFloat            lastMaxVentricularPressure_ = 0;
    TFloat            r1_                         = 0;
    TFloat            r2_                         = 0;
    TFloat            c_                            = 0;
    TFloat            c2_                           = 0;
    TFloat            b_                            = 0;
    TFloat            lastB_                        = 0;
    TFloat            d_                            = 0;
    TFloat            loadingPressure_              = 0;
    TFloat            ventricularResidualPressure_  = 0;
//    TFloat            c4_                           = 0;
//    TFloat            b4_                           = 0;
//    TFloat            d4_                           = 0;
    
    bool              minSystolicPressureSucceeded_ = false;
    bool              stepBack_                     = false;
    bool              shallPrepare_                 = true;
    
    TInt               relaxCnt_                     = 0;
    TInt               lastRelaxCnt_                 = 0;
    TInt               initMethod_                   = 1;
    TInt              minRelaxCnt_                = 0;
    TInt               lastValveState_             = 0;
    TInt               valveState_                 = 0;
    std::vector<TInt> preloadingMaterialIndices_    = {};
    
    typedef CBCirculatorySystemModel   Base;
};
#endif
