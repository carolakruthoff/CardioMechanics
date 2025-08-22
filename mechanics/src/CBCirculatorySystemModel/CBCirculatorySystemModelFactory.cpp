/*
 *  CBCirculatorySystemModelFactory.cpp
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 14.09.11.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#include "CBCirculatorySystemModelFactory.h"

#include "CBWindkessel3.h"
#include "CBExtension.h"

CBCirculatorySystemModel* CBCirculatorySystemModelFactory::New(CBCavity* cavity, ParameterMap* parameters)
{
    std::string type = parameters->Get<std::string>("Plugins.CirculatorySystem.Cavities.Cavity_" + std::to_string(cavity->GetIndex()) + ".Type","Windkessel3");
    CBCirculatorySystemModel* model=0;
    if(type == "Windkessel3")
        model = new CBWindkessel3(cavity, parameters);
    if(!model)
        model = CBExtensionLoadCirculatorySystemModel(cavity,parameters);
    if(!model)
        throw std::runtime_error("CBCirculatorySystemModel* CBCirculatorySystemModelFactory::New(CBCavity* cavity,ParameterMap* parameters): Model type " + type + " is unkown ");
    return model;
}



