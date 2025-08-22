/*
 *  CBPressureVolumeData.h
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 12.05.11.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#ifndef CB_PRESSURE_VOLUME_DATA
#define CB_PRESSURE_VOLUME_DATA

#include "Matrix3.h"

#include "CBElementCavity.h"
#include "CBCirculatorySystemModel.h"
#include "ParameterMap.h"

#include <string>
#include <vector>

/*
 * See:
 *
 * Total arterial inertance as the fourth element of the
 * windkessel model
 * Nikos Stergiopulos, Berend E. Westerhof and Nico Westerhof
 * Am J Physiol Heart Circ Physiol 276:H81-H88, 1999.
 */


using namespace math_pack;


class CBPressureVolumeData : public CBCirculatorySystemModel
{
public:
    CBPressureVolumeData(CBCavity* cavity, ParameterMap* parameters);
    TFloat Update(TFloat time, TFloat vol) override;
    TFloat GetPressure() override {return(ventricularPressure_); }
    TFloat CalcPressure(TFloat vol) override;
    TFloat CalcPressureDerivative(TFloat vol) override;
    void Prepare() override;
    double GetPreparationProgress() override {return time_/(preloadTime1_+preloadTime2_);}
    void StepBack() override;
    void WriteToFile(TFloat time) override;
    void InitInitialVolume() override;
protected:
private:
    void GetCorrespondingPressure(TFloat volume);
    TFloat ventricularPressure_;
    TFloat isovolumetricRelaxationCavityVolume_;
    TFloat atrialPressure_;
    TFloat aorticPressure_;
    TFloat volume_;
    TFloat initialVolume_;
    int    valveState_;
    TFloat initialAorticPressure_;
    TFloat lastVentricularPressure_;
    TFloat lastAorticPressure_;
    TFloat lastVolume_;
    int    lastValveState_;
    TFloat preloadTime1_;
    TFloat preloadTime2_;
    TFloat prevVolume_;
    TFloat lastPrevVolume_;
    TFloat bulkModulus_;
    TFloat dt_;
    TFloat prevDt_;
    TFloat lastPrevDt_;
    TInt minRelaxCnt_;
    TFloat minSystolicPressure_;
    TFloat maxVentricularPressure_;
    TFloat lastMaxVentricularPressure_;
    TFloat r1_;
    TFloat r2_;
    TFloat c_;
    TFloat c2_;
    TFloat b_;
    TFloat lastB_;
    TFloat d_;
    TFloat loadingPressure_;
    TFloat ventricularResidualPressure_;
    std::vector<TInt> preloadingMaterialIndices_;
    bool  minSystolicPressureSucceeded_ = false;
    bool  stepBack_ = false;
    bool shallPrepare_ = true;
    int   relaxCnt_;
    int   lastRelaxCnt_;
    int   initMethod_ = 1;
    std::vector<std::pair<double,double>> pv_;
    
    typedef CBCirculatorySystemModel   Base;
};

#endif

