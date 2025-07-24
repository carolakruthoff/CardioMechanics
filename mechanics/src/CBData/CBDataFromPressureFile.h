/*
 *  CBDataFromPressureFile.h
 *  CardioMechanics
 *
 *  Created by ek717 on 25.04.2020
 *  Copyright 2010 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#ifndef CB_DATA_FROM_PRESSURE_FILE_H
#define CB_DATA_FROM_PRESSURE_FILE_H

#include <stdint.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

#include "ParameterMap.h"

#include "CBData.h"
#include "DCType.h"

class CBData;


class CBDataFromPressureFile  : public CBData
{
public:
    CBDataFromPressureFile();
    CBDataFromPressureFile(ParameterMap* parameters, std::string parameterKey);
    void Init(std::string filename);
    virtual TFloat Get(TFloat time, TInt index);

protected:

private:
    //virtual void LoadFileList(std::string filename);
    bool LoadDataSet(std::string filename); //TFloat time);

    std::vector<std::string> SplitString(std::string input);

    std::vector<TFloat>                         dataTimeSteps_;
    std::vector<TFloat>                         dataPressureAtTimeStep_;
};


#endif
