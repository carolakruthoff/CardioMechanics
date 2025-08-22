/*
 *  CBCavity.h
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 31.03.11.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#ifndef CB_CAVITY_H
#define CB_CAVITY_H

#include "Matrix3.h"

#include "ParameterMap.h"
#include "CBCirculatorySystemModel.h"
#include "CBCirculatorySystemModelFactory.h"
#include "CBElementCavity.h"


class CBCirculatorySystemModel;


using namespace math_pack;

// TODO: What is the difference between CBCirculationCavity and CBCavity?

class CBCavity
{
public:
  CBCavity(){}
  ~CBCavity();
  void InsertNextElement(CBElementCavity* element);
  void SetIndex(TInt index);
  TInt GetIndex();
  void SetInitialVolume(TFloat initialVolume);
  TFloat GetInitialVolume();
  void SetVolume(TFloat vol);
  TFloat GetVolume();
  void Update(TFloat time, TFloat vol);
  void InitPressure();
  void SetAdapter(CBElementAdapter* adapter){adapter_ = adapter;}
  CBElementAdapter* GetAdapter(){return adapter_;}
  bool IsActive(){return model_->IsActive();}
  CBStatus GetStatus();
  void Prepare();
  double GetPreparationProgress(){return model_->GetPreparationProgress();}
  TFloat GetPressure();
  TFloat CalcPressureDerivative(TFloat vol);
  TFloat CalcPressure(TFloat vol);
  void StepBack();
  void WriteToFile(TFloat time);
  void InitCirculatorySystemModel(ParameterMap* parameters);
  void SetCirculatorySystemModel(CBCirculatorySystemModel* model);
  CBCirculatorySystemModel* GetCirculatorySystemModel();
  const std::vector<CBElementCavity*> GetCavityElements() {return elements_;}
  
protected:
private:
  std::vector<CBElementCavity* > elements_;
  TInt                           index_ = 0;
  TFloat                         initialVolume_=0;
  TFloat                         volume_ = 0;
  TFloat                         pressure_ = 0;
  CBElementAdapter*              adapter_ = 0;

  CBCirculatorySystemModel*      model_ = 0;
  bool                           isIndexSet_=false;
};

#endif
