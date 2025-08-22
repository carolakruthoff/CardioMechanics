/*
 *  CBCirculatorySystemModelFactory.h
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 16.06.11.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#ifndef CB_CIRCULATORY_SYSTEM_MODEL_FACTORY_H
#define CB_CIRCULATORY_SYSTEM_MODEL_FACTORY_H

#include "Matrix3.h"

#include "ParameterMap.h"

#include "CBCavity.h"


class CBCirculatorySystemModelFactory
{
public:
    static CBCirculatorySystemModel* New(CBCavity* cavity, ParameterMap* parameters);
protected:
private:
    CBCirculatorySystemModelFactory();
};

#endif