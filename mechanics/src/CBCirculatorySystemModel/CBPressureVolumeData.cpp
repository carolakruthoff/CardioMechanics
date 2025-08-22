/*
 *  CBPressureVolumeData.cpp
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 06.04.13.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#include <cmath>

#include "CBSolver.h"
#include "CBPressureVolumeData.h"
#include "CBCavity.h"

void CBPressureVolumeData::InitInitialVolume()
{
    initialVolume_ = Base::cavity_->GetInitialVolume();
    volume_        = initialVolume_;
}


void CBPressureVolumeData::WriteToFile(TFloat time)
{
    if(Base::filename_ != "")
    {
        Base::file_.open(Base::filename_.c_str(), std::ios::app);
        
        if(!file_.good())
            throw std::runtime_error("void CBPressureVolumeData::WriteToFile(TFloat time): Couldn't create " + Base::filename_ + ".");
        
        Base::file_ << time << "\t\t" << volume_ << "\t\t" << ventricularPressure_ << "\t\t\t\t" << aorticPressure_ << std::endl;
        Base::file_.close();
    }
}


void CBPressureVolumeData::Prepare()
{
    if(Base::isActive_ && shallPrepare_)
    {
        valveState_ = lastValveState_ = -2;
        status_     = CBStatus::PREPARING_SIMULATION;
    }
    else
    {
        status_     = CBStatus::DACCORD;
        valveState_ = -1;
    }
}


void CBPressureVolumeData::StepBack()
{
    // Step back in time if solution was not accepted
    
    stepBack_ = true;
    
    ventricularPressure_ = lastVentricularPressure_;
    aorticPressure_      = lastAorticPressure_;
    volume_              = lastVolume_;
    prevVolume_          = lastPrevVolume_;
    
    if(lastValveState_ != 1 && lastValveState_ != -4 && lastValveState_ != -2)
        valveState_ = lastValveState_;
    
    relaxCnt_ = lastRelaxCnt_;
    maxVentricularPressure_ = lastMaxVentricularPressure_;
    Base::time_            -= dt_;
    prevDt_                 = lastPrevDt_;
    b_ = lastB_;
}


TFloat CBPressureVolumeData::Update(TFloat time, TFloat currentVolume)
{
    if(Base::isActive_ && !(status_ == CBStatus::WAITING))
    {
        if(!stepBack_)
        {
            lastB_      = b_;
            lastVolume_ = volume_;
            lastVentricularPressure_    = ventricularPressure_;
            lastAorticPressure_         = aorticPressure_;
            lastValveState_             = valveState_;
            lastRelaxCnt_               = relaxCnt_;
            lastMaxVentricularPressure_ = maxVentricularPressure_;
            lastPrevVolume_             = prevVolume_;
            lastPrevDt_ = prevDt_;
        }
        
        stepBack_ = false;
        dt_ = cavity_->GetAdapter()->GetSolver()->GetTiming().GetTimeStep();
        
        if(prevDt_ == 0)
            prevDt_ = dt_;
        
        ventricularPressure_ = CalcPressure(currentVolume);
        
        if(status_ == CBStatus::PREPARING_SIMULATION)
        {
            time = Base::time_ + dt_;
            TFloat db = (bulkModulus_ / (preloadTime2_))*2;
            
            if(valveState_ == -2)
            {
                if(time > preloadTime1_/2.0)
                    loadingPressure_ = 0;
                else
                {
                    if(preloadingMaterialIndices_.size() == 0)
                        cavity_->GetAdapter()->GetSolver()->RelaxElementsAndBases();
                    else
                        for(auto i : preloadingMaterialIndices_)
                            cavity_->GetAdapter()->GetSolver()->RelaxElementsAndBasesByMaterialIndex(i);
                }
                
                if(time > preloadTime1_)
                {
                    if(initMethod_ == 1)
                        valveState_ = -1;
                    else
                        if(initMethod_ == 2)
                            valveState_ = -4;
                        else
                            throw std::runtime_error("TFloat CBPressureVolumeData::Update(TFloat time, TFloat currentVolume): Unkown initialization method");
                    time        = 0;
                }
            }
            else
            {
                switch (initMethod_)
                {
                    case 1:
                        
                        if(b_ < bulkModulus_)
                            b_ = db*time;
                        else
                            b_ = bulkModulus_;
                        if(time > preloadTime2_)
                        {
                            status_ = CBStatus::DACCORD;
                            ventricularResidualPressure_ = ventricularPressure_;
                            time = 0;
                        }
                        break;
                    case 2:
                    {
                        
                        TFloat dv =2*(initialAorticPressure_ / preloadTime2_);
                        
                        if(ventricularPressure_ < initialAorticPressure_)
                            ventricularPressure_ = dv * time;
                        else
                            ventricularPressure_ = initialAorticPressure_;
                        
                        if(time > preloadTime2_)
                        {
                            valveState_ = 2;
                            time = 0;
                            status_ = CBStatus::DACCORD;
                        }
                        break;
                    }
                    default:
                        throw std::runtime_error("TFloat CBPressureVolumeData::Update(TFloat time, TFloat currentVolume): Something went really wrong !!!");
                        break;
                }
            }
        }
        else
        {
            if(valveState_ == -1)
                valveState_ = 1;
            
            if(ventricularPressure_ > maxVentricularPressure_)
            {
                maxVentricularPressure_ = ventricularPressure_;
                relaxCnt_ = 0;
            }
            else
                if(valveState_ == 2)
                    relaxCnt_++;
            
            if((valveState_ == 0) || (valveState_ == 1) || (valveState_ == 3))
                aorticPressure_ -= (aorticPressure_ * dt_) / (r2_ * c_);
            else if(valveState_ == 2)
                aorticPressure_ = ventricularPressure_;
            
            if((valveState_ == 1) && (ventricularPressure_ > aorticPressure_))
                valveState_ = 2;
            
            if((valveState_ == 2) && (ventricularPressure_ < maxVentricularPressure_) && (relaxCnt_ > minRelaxCnt_) && minSystolicPressureSucceeded_ && (currentVolume > volume_))
            {
                valveState_ = 3;
                isovolumetricRelaxationCavityVolume_ = currentVolume;
            }
            
            if((valveState_ == 3) && (ventricularPressure_ < atrialPressure_))
                valveState_ = 4;
            
            if((valveState_ == 4) && (ventricularPressure_ <= ventricularResidualPressure_))
            {
                valveState_          = 0;
                ventricularPressure_ = ventricularResidualPressure_;
            }
            
            if(ventricularPressure_ > minSystolicPressure_)
                minSystolicPressureSucceeded_ = true;
        }
        
        prevVolume_ = volume_;
        prevDt_     = dt_;
        volume_     = currentVolume;
        Base::time_ = time;
        
        DCCtrl::debug << "\n-----------------------------------\nCavity: " <<  Base::cavity_->GetIndex() << "\nStage: ";
        if(status_==CBStatus::PREPARING_SIMULATION)
            DCCtrl::debug << std::string("Initializing [")+std::to_string(initMethod_)+std::string("]");
        else
            DCCtrl::debug << std::string("Running");
        DCCtrl::debug << "\nValve state: "<< valveState_ << "\nTime: " << time << "s\nBulk modulus: " << b_ << "\nCavity Pressure: " << ventricularPressure_/133.322 << " mmHg\nAfter Load Pressure: " << aorticPressure_/133.322 << " mmHg\nVolume: " << volume_ * 1.0e6 << " ml\nInitial volume: " << initialVolume_ * 1.0e6 << " ml" << "\n-----------------------------------\n";
        
        return(ventricularPressure_);
    }
    else
    {
        volume_ = currentVolume;
        return(0);
    }
}


TFloat CBPressureVolumeData::CalcPressure(TFloat currentVolume)
{
    if(Base::isActive_ && !(status_ == CBStatus::WAITING))
        switch(valveState_)
    {
        case -4:
            return ventricularPressure_;
            break;
        case -3:
            return(0);
            
            break;
            
        case -2:
            return(-loadingPressure_);
            
            break;
            
        case -1:
            return(-b_ * (currentVolume / initialVolume_ - 1.0));
            
            break;
            
        case 0:
            return(ventricularResidualPressure_);
            
            break;
            
        case 1:
            return(-b_ * (currentVolume / initialVolume_ - 1.0));
            
            break;
            
        case 2:
            return(aorticPressure_ - (aorticPressure_ * dt_) / (r2_ * c_) - (1 + r1_/r2_)*(currentVolume - volume_) / c_ - r1_ * ((currentVolume - volume_) / dt_ - (volume_ - prevVolume_) / prevDt_));
            break;
            
        case 3:
            return(-bulkModulus_ * (currentVolume / isovolumetricRelaxationCavityVolume_ - 1.0));
            break;
            
        case 4:
            return(ventricularPressure_ - (d_*(ventricularPressure_-ventricularResidualPressure_) * dt_) - (currentVolume -volume_) / c2_);
            
            break;
            
        default:
            throw std::runtime_error("void CBPressureVolumeData::CalcPressure(TFloat currentVolume): unsupported valve state");
    }
    
    
    else
        return(0);
}


TFloat CBPressureVolumeData::CalcPressureDerivative(TFloat vol)
{
    if(Base::isActive_ && !(status_ == CBStatus::WAITING))
        switch(valveState_)
    {
        case -4:
            return(0);
            
            break;
            
        case -3:
            return(0);
            
            break;
            
        case -2:
            return(0);
            
            break;
            
        case -1:
            return(-b_  / initialVolume_);
            
            break;
            
        case 0:
            return(0);
            
            break;
            
        case 1:
            return(-b_  / initialVolume_);
            
            break;
            
        case 2:
            return(-(1+r1_/r2_) / c_ - r1_/dt_);
            
            break;
            
        case 3:
            return(-bulkModulus_  / isovolumetricRelaxationCavityVolume_);
            
            break;
            
        case 4:
            return(-1 / c2_);
            
            break;
            
        default:
            throw std::runtime_error("void CBPressureVolumeData::CalcPressure(TFloat currentVolume): unsupported valve state");
    }
    
    
    else
        return(0);
}



void CBPressureVolumeData::GetCorrespondingPressure(TFloat volume)
{
    
}

CBPressureVolumeData::CBPressureVolumeData(CBCavity* cavity, ParameterMap* parameters) : CBCirculatorySystemModel(cavity, parameters)
{
    std::string index = std::to_string(Base::cavity_->GetIndex());
    
    
    std::string filename = parameters->Get<std::string>("Plugins.CirculatorySystem.Cavities.Cavity_" + index + ".File","");
    
    
    Base::isActive_ = parameters->Get<bool>("Plugins.CirculatorySystem.Cavities.Cavity_" + index + ".Active", false);
    
    if(Base::isActive_ && (Base::filename_ != ""))
    {
        std::ifstream file(filename.c_str());
        if(!file.good())
            throw std::runtime_error("void void CBPressureVolumeData::CBPressureVolumeData(): File " + filename + " does not exists");
        std::string str;
        getline(file,str); // ignore first line
        
        while (getline(file,str))
        {
            std::stringstream ss(str);
            TFloat time, vol, ventricularPressure;
            ss >> time >> vol >> ventricularPressure;
            pv_.push_back({vol,ventricularPressure});
        }
        
        file.close();
    }
    else
    {
        aorticPressure_      = 0;
        ventricularPressure_ = 0;
    }
    
    atrialPressure_ = parameters->Get<TFloat>("Plugins.CirculatorySystem.Cavities.Cavity_" + index + ".ClosingPressure", 2000);
    ventricularResidualPressure_ = parameters->Get<TFloat>("Plugins.CirculatorySystem.Cavities.Cavity_" + index + ".VentricularResidualPressure", 750);
    minSystolicPressure_         = parameters->Get<TFloat>("Plugins.CirculatorySystem.Cavities.Cavity_" + index + ".MinSystolicPressure", 10000.0);
    
    bulkModulus_ = parameters->Get<TFloat>("Plugins.CirculatorySystem.Cavities.Cavity_" + index + ".PressureVolumeData.BulkModulus", 1e7);
    c_           = parameters->Get<TFloat>("Plugins.CirculatorySystem.Cavities.Cavity_" + index + ".PressureVolumeData.C", 1e-8);
    c2_          = parameters->Get<TFloat>("Plugins.CirculatorySystem.Cavities.Cavity_" + index + ".PressureVolumeData.C2", 1e-7);
    r1_          = parameters->Get<TFloat>("Plugins.CirculatorySystem.Cavities.Cavity_" + index + ".PressureVolumeData.R1", 1.4e8);
    r2_          = parameters->Get<TFloat>("Plugins.CirculatorySystem.Cavities.Cavity_" + index + ".PressureVolumeData.R2", 1.4e8);
    d_           = parameters->Get<TFloat>("Plugins.CirculatorySystem.Cavities.Cavity_" + index + ".PressureVolumeData.D", 12);
    
    minRelaxCnt_     = parameters->Get<TFloat>("Plugins.CirculatorySystem.Cavities.Cavity_" + index + ".MinRelaxCnt", 5);
    
    preloadingMaterialIndices_ = parameters->GetArray<TInt>("Plugins.CirculatorySystem.Cavities.Cavity_" + index + ".PreloadMaterialIndices",{});
    if(Base::isActive_)
    {
        aorticPressure_      = parameters->Get<TFloat>("Plugins.CirculatorySystem.Cavities.Cavity_" + index + ".AorticPressure", 9000);
        ventricularPressure_ = parameters->Get<TFloat>("Plugins.CirculatorySystem.Cavities.Cavity_" + index + ".VentricularPressure", 0);
    }
    else
    {
        aorticPressure_      = 0;
        ventricularPressure_ = 0;
    }
    
    lastVentricularPressure_ = ventricularPressure_;
    
    initialAorticPressure_     = aorticPressure_;
    lastAorticPressure_      = aorticPressure_;
    
    valveState_ =  parameters->Get<TFloat>("Plugins.CirculatorySystem.Cavities.Cavity_" + index + ".ValveState", -3);
    shallPrepare_ = parameters->Get<bool>("Plugins.CirculatorySystem.Cavities.Cavity_" + index + ".Initialize", true);
    loadingPressure_ = parameters->Get<TFloat>("Plugins.CirculatorySystem.Cavities.Cavity_" + index + ".LoadingPressure", 600);
    initMethod_ = parameters->Get<TFloat>("Plugins.CirculatorySystem.Cavities.Cavity_" + index + ".InitMethod", 1);
    
    preloadTime1_    = parameters->Get<TFloat>("Plugins.CirculatorySystem.PreloadTime1", 0.02);
    preloadTime2_    = parameters->Get<TFloat>("Plugins.CirculatorySystem.PreloadTime2", 0.1);
    
    lastValveState_ = valveState_;
    
    b_     = 1;
    lastB_ = 1;
    
    relaxCnt_     = 0;
    lastRelaxCnt_ = 0;
    
    volume_ = 0;
    lastVolume_ = volume_;
    prevVolume_ = volume_;
    lastPrevVolume_ = volume_;
    
    dt_ = 0;
    
    minSystolicPressureSucceeded_ = false;
    
    if(Base::filename_ != "")
    {
        Base::file_.open(Base::filename_.c_str());
        Base::file_ << "Time \t\t Volume \t\t  Ventricular Pressure \t\t Aortic Pressure \n";
        Base::file_.close();
    }
}
