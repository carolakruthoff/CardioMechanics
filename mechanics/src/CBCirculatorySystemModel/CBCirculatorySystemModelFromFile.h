//
//  CBCirculatorySystemModelFromFile.h
//  CardioMechanics_unstable
//
//  Created by Thomas Fritz on 07.03.13.
//
//



#ifndef CB_CIRCULATORY_SYSTEM_MODEL_FROM_FILE_H
#define CB_CIRCULATORY_SYSTEM_MODEL_FROM_FILE_H

#include "CBCirculatorySystemModelFromFile.h"
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


class CBCirculatorySystemModelFromFile : public CBCirculatorySystemModel
{
public:
    CBCirculatorySystemModelFromFile(CBCavity* cavity, ParameterMap* parameters);
    TFloat Update(TFloat time, TFloat vol);
    TFloat GetPressure(){return(ventricularPressure_); }
    TFloat CalcPressure(TFloat vol);
    TFloat CalcPressureDerivative(TFloat vol);
    void WriteToFile(TFloat time);
    void InitInitialVolume();
    
protected:
private:
    
    std::vector<TFloat> pressureAtTimeStep_;
    std::vector<TFloat> timeSteps_;
    TInt currentTimeStep_;
    TFloat ventricularPressure_;
    TFloat aorticPressure_;
    TFloat volume_;
    TFloat initialVolume_;
    
    typedef CBCirculatorySystemModel   Base;
};

#endif
