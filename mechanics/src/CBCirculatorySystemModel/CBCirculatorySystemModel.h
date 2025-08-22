/*
 *  CBCirculatorySystemModel.h
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 12.05.11.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#ifndef CB_CIRCULATORY_SYSTEM_MODEL_H
#define CB_CIRCULATORY_SYSTEM_MODEL_H

#include "Matrix3.h"

#include "CBStatus.h"
#include "DCType.h"

#include "CBElementCavity.h"
#include "ParameterMap.h"

using namespace math_pack;

class CBCavity;

class CBCirculatorySystemModel
{
public:
    CBCirculatorySystemModel(CBCavity* cavity, ParameterMap* parameters);
    virtual ~CBCirculatorySystemModel()
    {
        file_.close();
    }
    virtual void Init(ParameterMap* parameters){}
    virtual TFloat Update(TFloat time, TFloat vol) = 0;
    virtual TFloat GetPressure() = 0;
    virtual TFloat CalcPressure(TFloat vol) = 0;
    virtual TFloat CalcPressureDerivative(TFloat vol) = 0;
    virtual void WriteToFile(TFloat time) = 0;
    virtual void InitInitialVolume() = 0;
    virtual void Prepare() { std::runtime_error("Error: Function CBCirculatorySystemModel::Prepare() is not implemented."); }
    virtual double GetPreparationProgress(){return -1.0;}
    virtual void StepBack() { std::runtime_error("Error: Function CBCirculatorySystemModel::StepBack() is not implemented."); }
    bool IsActive(){return(isActive_);}
    CBStatus GetStatus(){return(status_);}

protected:
    bool          isActive_;
    CBCavity*     cavity_;
    std::string   filename_;
    std::ofstream file_;
    CBStatus      status_ = CBStatus::WAITING;
    TFloat        time_;

private:
};
#endif
