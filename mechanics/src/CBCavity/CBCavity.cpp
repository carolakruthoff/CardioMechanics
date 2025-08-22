/*
 *  CBCavity.cpp
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 14.09.11.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#include "CBCavity.h"
//myTest...


CBCavity::~CBCavity()
{
  for( auto& it: elements_)
        delete it;
    if(model_)
        delete model_;
    model_ = 0;
}

void CBCavity::InsertNextElement(CBElementCavity* element)
{
    elements_.push_back(element); 
}



void CBCavity::SetIndex(TInt index)
{
    index_      = index;
    isIndexSet_ = true;
}


TInt CBCavity::GetIndex()
{
    return(index_); 
}


void CBCavity::SetInitialVolume(TFloat initialVolume)
{
    initialVolume_ = initialVolume;
    model_->InitInitialVolume();
}


TFloat CBCavity::GetInitialVolume()
{
    return(initialVolume_); 
}


void CBCavity::SetVolume(TFloat vol)
{
    volume_ = vol;
}


TFloat CBCavity::GetVolume()
{
    return(volume_);
}


void CBCavity::Update(TFloat time, TFloat vol)
{
    volume_   = vol;
    pressure_ = model_->Update(time, vol);
}


void CBCavity::InitPressure()
{
    pressure_ = model_->GetPressure();
}

void CBCavity::Prepare()
{
    model_->Prepare();
}

void CBCavity::StepBack()
{
    model_->StepBack();
}

CBStatus CBCavity::GetStatus()
{
    return(model_->GetStatus());
}


TFloat CBCavity::GetPressure()
{
    return(pressure_);
}


TFloat CBCavity::CalcPressureDerivative(TFloat vol)
{
    return(model_->CalcPressureDerivative(vol));
}


TFloat CBCavity::CalcPressure(TFloat vol)
{
    return(model_->CalcPressure(vol));
}



void CBCavity::WriteToFile(TFloat time)
{
    model_->WriteToFile(time);
}


void CBCavity::SetCirculatorySystemModel(CBCirculatorySystemModel* model)
{
    model_ = model;
}


CBCirculatorySystemModel* CBCavity::GetCirculatorySystemModel()
{
    return(model_);
}



void CBCavity::InitCirculatorySystemModel(ParameterMap* parameters)
{
    if(isIndexSet_)
        model_ = CBCirculatorySystemModelFactory::New(this, parameters);
    else
        throw std::runtime_error("class CBCavity::InitCirculatorySystemModel(ParameterMap* parameters): SetIndex(TInt index) has to be run first");
}

