//
//  CBCirculatorySystemModelFromFile.cpp
//  CardioMechanics_unstable
//
//  Created by Thomas Fritz on 07.03.13.
//
//

#include "CBCirculatorySystemModelFromFile.h"
#include "CBCavity.h"

void CBCirculatorySystemModelFromFile::InitInitialVolume()
{
    initialVolume_ = Base::cavity_->GetInitialVolume();
    volume_        = initialVolume_;
}


void CBCirculatorySystemModelFromFile::WriteToFile(TFloat time)
{
    // Do nothing
}


TFloat CBCirculatorySystemModelFromFile::Update(TFloat time, TFloat currentVolume)
{
    if(Base::isActive_)
    {
        if(time > timeSteps_.at(0))
        {
            while(time > timeSteps_.at(currentTimeStep_+1) && currentTimeStep_ < timeSteps_.size()-2)
                currentTimeStep_++;

            TFloat t1 = timeSteps_.at(currentTimeStep_);
            TFloat p1 = pressureAtTimeStep_.at(currentTimeStep_);
            TFloat t2,p2;
            if(currentTimeStep_ < timeSteps_.size()-1)
            {
                t2 = timeSteps_.at(currentTimeStep_+1);
                p2 = pressureAtTimeStep_.at(currentTimeStep_+1);
            }
            else
            {
                t2 = -timeSteps_.at(currentTimeStep_-1);
                p2 = -pressureAtTimeStep_.at(currentTimeStep_-1);
            }
            ventricularPressure_ = p1 + (p2 - p1)/(t2-t1)* (time - t1);
        }
        else
            ventricularPressure_ = 0;
    }
    else
        ventricularPressure_ = 0;
    std::cout << "FromFile" << time << " " << currentVolume << " " << ventricularPressure_ << "\n";
    return(ventricularPressure_);
}


TFloat CBCirculatorySystemModelFromFile::CalcPressure(TFloat currentVolume)
{
    if(Base::isActive_)
        return(ventricularPressure_);
    else
        return(0);
}


TFloat CBCirculatorySystemModelFromFile::CalcPressureDerivative(TFloat vol)
{
    return(0);
}


CBCirculatorySystemModelFromFile::CBCirculatorySystemModelFromFile(CBCavity* cavity, ParameterMap* parameters) : CBCirculatorySystemModel(cavity, parameters)
{
    std::string index = std::to_string(Base::cavity_->GetIndex());
    Base::isActive_ = parameters->Get<bool>("Plugins.CirculatorySystem.Cavities.Cavity_" + index + ".Active", true);
    std::string filename = parameters->Get<std::string>("Plugins.CirculatorySystem.Cavities.Cavity_" + index + ".FromFile.Filename","");

    if(Base::isActive_ && (Base::filename_ != ""))
    {
        std::ifstream file(filename.c_str());
        if(!file.good())
            throw std::runtime_error("CBCirculatorySystemModelFromFile::CBCirculatorySystemModelFromFile(CBCavity* cavity, ParameterMap* parameters, TFloat startTime) : CBCirculatorySystemModel(cavity, parameters, startTime): File " + filename + " does not exists");
        currentTimeStep_ = 0;
        std::string str;
        getline(file,str); // ignore first line
        while (getline(file,str))
        {
            std::stringstream ss(str);
            TFloat time, vol, ventricularPressure;
            ss >> time >> vol >> ventricularPressure;
            timeSteps_.push_back(time);
            pressureAtTimeStep_.push_back(ventricularPressure);
        }
        file.close();
    }
    else
    {
        aorticPressure_      = 0;
        ventricularPressure_ = 0;
    }
}