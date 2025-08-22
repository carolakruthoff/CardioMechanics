/*
 *  CBCirculatorySystemModel.cpp
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 14.09.11.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */


#include "filesystem.h"

#include "CBCirculatorySystemModel.h"
#include "CBCavity.h"



CBCirculatorySystemModel::CBCirculatorySystemModel(CBCavity* cavity, ParameterMap* parameters) : time_(0)
{
    cavity_   = cavity;
    std::string index = std::to_string(cavity_->GetIndex());
    filename_ = parameters->Get<std::string>("Plugins.CirculatorySystem.Cavities.Cavity_" + index + ".ExportFile", "");
    
    if(filename_!="")
    {

        size_t pos = filename_.find_last_of("/\\");
        if(pos != std::string::npos)
        {
            std::string dir        = filename_.substr(0, pos);
            frizzle::filesystem::CreateDirectory(dir);
        }
    }
}

