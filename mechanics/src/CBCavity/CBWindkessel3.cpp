/*
 *  CBWindkessel3.cpp
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 06.04.13.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#include <cmath>

#include "CBSolver.h"
#include "CBWindkessel3.h"
#include "CBCavity.h"

void CBWindkessel3::InitInitialVolume()
{
    initialVolume_ = Base::cavity_->GetInitialVolume();
    volume_        = initialVolume_;
}


void CBWindkessel3::WriteToFile(TFloat time)
{
    if(Base::filename_ != "")
    {
        Base::file_.open(Base::filename_.c_str(), std::ios::app);

        if(!file_.good())
            throw std::runtime_error("void CBWindkessel3::WriteToFile(TFloat time): Couldn't create " + Base::filename_ + ".");

        Base::file_ << time << "\t\t" << volume_ << "\t\t" << ventricularPressure_ << "\t\t\t\t" << aorticPressure_ << "\t\t\t\t" << (volume_ - lastVolume_)/dt_  << "\t\t\t\t" << (ventricularPressure_- lastVentricularPressure_)/dt_ << "\t\t\t\t" << (lastVentricularPressure_- lastAorticPressure_)/dt_ <<  "\t\t\t\t" << volumePressureWork_ << "\t\t" << valveState_ << std::endl;
        Base::file_.close();
    }
}


void CBWindkessel3::Prepare()
{
    if(Base::isActive_ && status_ == CBStatus::WAITING)
    {
        if(shallPrepare_)
        {
            valveState_ = -2;
            lastValveState_ = -2;
            status_     = CBStatus::PREPARING_SIMULATION;
        }
        else
        {
            b_ = bulkModulus_;
            valveState_ = -1;
            status_     = CBStatus::DACCORD;
        }
    }
}


void CBWindkessel3::StepBack()
{
    // Step back in time if solution was not accepted
    if(Base::IsActive() && status_ != CBStatus::WAITING)
    {
        stepBack_ = true;

        ventricularPressure_ = lastVentricularPressure_;
        volumePressureWork_  = lastVolumePressureWork_;
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
}


TFloat CBWindkessel3::Update(TFloat time, TFloat currentVolume)
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
            lastVolumePressureWork_ = volumePressureWork_;
        }

        stepBack_ = false;
        dt_ = cavity_->GetAdapter()->GetSolver()->GetTiming().GetTimeStep();

        if(prevDt_ == 0)
            prevDt_ = dt_;

        // TODO: call CalcPressure AFTER valve states are correctly set!
        ventricularPressure_ = CalcPressure(currentVolume);

        if(status_ == CBStatus::PREPARING_SIMULATION)
        {
            time = Base::time_ + dt_;
            if(valveState_ == -2)
            {
                if(time > 0.5*preloadTime1_)
                    loadingPressure_ = 0;
                else
                {
                    if(preloadingMaterialIndices_.size() == 0)
                        cavity_->GetAdapter()->GetSolver()->RelaxElementsAndBases();
                    else
                        for(auto i : preloadingMaterialIndices_) {
                            // set zero stress to the current state
                            cavity_->GetAdapter()->GetSolver()->RelaxElementsAndBasesByMaterialIndex(i);
                        }
                }

                if(time > preloadTime1_)
                {
                    if(initMethod_ == 1)
                        valveState_ = -1;
                    else
                        if(initMethod_ == 2)
                            valveState_ = -4;
                        else
                            throw std::runtime_error("TFloat CBWindkessel3::Update(TFloat time, TFloat currentVolume): Unkown initialization method");
                }
            }
            else
            {
                switch (initMethod_)
                {
                    // TODO: How is pressure change achieved / computed in initMethod_=1 ?
                    case 1:
                        {
                            if (time < preloadTime1_) {
                                b_ = 1;
                            }
                            else if (time < preloadTime1_ + 0.5*preloadTime2_) {
                                TFloat t = (time - preloadTime1_) / (0.5*preloadTime2_);
                                b_ = 1 + (bulkModulus_-1)*(10+(6*t-15)*t)*pow(t, 3.);
                            }
                            else {
                                b_ = bulkModulus_;
                            }

                            if(time > preloadTime1_ + preloadTime2_)
                            {
                                time = 0;
                                valveState_ = 1;
                                ventricularResidualPressure_ = ventricularPressure_;
                                status_ = CBStatus::DACCORD;
                            }
                            break;
                        }
                    case 2:
                        {
                            if ( preloadTime1_<=time && time<preloadTime1_+0.5*preloadTime2_ ) {
                                TFloat t = (time - preloadTime1_) / (0.5*preloadTime2_);
                                ventricularPressure_ = initialAorticPressure_*t;
                                ventricularPressure_ = initialAorticPressure_*(10+(6*t-15)*t)*pow(t,3.);
                            }
                            else
                                ventricularPressure_ = initialAorticPressure_;

                            if(time > preloadTime1_ + preloadTime2_)
                            {
                                valveState_ = 2;
                                time = 0;
                                status_ = CBStatus::DACCORD;
                            }
                            break;
                        }
                    default:
                        throw std::runtime_error("TFloat CBWindkessel3::Update(TFloat time, TFloat currentVolume): Something went really wrong !!!");
                        break;
                }
            }

        }
        else //  status != preparing_simulation
        {

            if ( (floor(time/beatInterval_) != floor(Base::time_/beatInterval_)) )
            {
                // WARNING: Replace this by DCCtrl::Print
                std::cout << "Beginning with a new heartbeat!" << std::endl;
                b_ = bulkModulus_;
                initialVolume_              = Base::cavity_->GetVolume();
                ventricularPressure_        = lastVentricularPressure_;
                initialAorticPressure_      = aorticPressure_;
                initialVentricularPressure_ = ventricularPressure_;
                valveState_             = 1;
                maxVentricularPressure_ = 0;
                relaxCnt_               = 0;
                volumePressureWork_     = 0;
            }

            // TODO: first define the final valve state, then compute the pressure (not mixed)!
            // (maybe use lastVentricularPressure_ ?)

            if(valveState_ == -1)
            {
                valveState_ = 1;
                volumePressureWork_ = 0;
            }
            if(valveState_ > 0 && valveState_ < 3 && currentVolume < lastVolume_)
            {
                volumePressureWork_ += (lastVolume_ - currentVolume)*ventricularPressure_;
            }
            if(valveState_ == 2 && ventricularPressure_ > maxVentricularPressure_)
            {
                maxVentricularPressure_ = ventricularPressure_;
                relaxCnt_ = 0;
            }
            else
                if(valveState_ == 2 || valveState_ == 3)
                    relaxCnt_++;

            if((valveState_ == 0) || (valveState_ == 1) || (valveState_ == 3) || (valveState_ == 4))
                aorticPressure_ -= (aorticPressure_ * dt_) / (r2_ * c_);
            else if(valveState_ == 2)
                aorticPressure_ = ventricularPressure_;

            if((valveState_ == 1) && (ventricularPressure_ > aorticPressure_)) {
                aorticPressure_ = ventricularPressure_;
                valveState_ = 2;
            }

            if((valveState_ == 2) && (ventricularPressure_ < maxVentricularPressure_) && (relaxCnt_ > minRelaxCnt_) && minSystolicPressureSucceeded_ && (currentVolume > volume_))
            {
                valveState_ = 3;
                isovolumetricRelaxationCavityVolume_ = currentVolume;
                relaxCnt_ = 2;
            }

            if((valveState_ == 3) && (ventricularPressure_ < atrialPressure_) && (relaxCnt_ > minRelaxCnt_))
            {
                valveState_ = 4;
                //                b4_ = ((ventricularPressure_ - lastVentricularPressure_)/(dt_*1000))/(ventricularResidualPressure_ - ventricularPressure_);
                //                c4_ = 1.6;
                //                d4_ = -(log(ventricularPressure_ - ventricularResidualPressure_)/b4_ + pow((time*1000),c4_));
            }

            if((valveState_ == 4) && (ventricularPressure_ <= ventricularResidualPressure_))
            {
                valveState_          = 0;
                ventricularPressure_ = ventricularResidualPressure_;
                minSystolicPressureSucceeded_ = false;
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
            DCCtrl::debug << std::string("Initializing [")+std::to_string(initMethod_)+std::string("]\n");
        else
            DCCtrl::debug << std::string("Running\n");

        DCCtrl::debug << "Valve state: "<< valveState_ << "\n"
            << "Time: " << time << "s\n"
            << "Bulk modulus: " << b_ << "\n"
            << "Cavity Pressure: " << ventricularPressure_/133.322 << "mmHg\n"
            << "Maximal Cavity Pressure: " << maxVentricularPressure_ /133.322 << "mmHg\n"
            << "Required Minimal Systolic Pressure: " << minSystolicPressure_/133.322 << "mmHg\n"
            << "Relaxation Counter: " << relaxCnt_<< "\n"
            << "After Load Pressure: " << aorticPressure_/133.322 << "mmHg\n"
            << "Volume: " << volume_ * 1.0e6 << " ml\n"
            << "Initial volume: " << initialVolume_ * 1.0e6 << " ml"
            << "\n-----------------------------------\n";

        return(ventricularPressure_);
    }
    else
    {
        volume_ = currentVolume;
        return(0);
    }
}


TFloat CBWindkessel3::CalcPressure(TFloat currentVolume)
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
                return -loadingPressure_;

                break;

            case -1:
                return( b_ * ((initialVolume_ - currentVolume) / initialVolume_) );

                break;

            case 0:
                return(ventricularResidualPressure_);

                break;

            case 1:
                return( b_ * ((initialVolume_ - currentVolume) / initialVolume_) );

                break;

            case 2:
                return(aorticPressure_ - (aorticPressure_ * dt_) / (r2_ * c_) - (1 + r1_/r2_)*(currentVolume - volume_) / c_ - r1_ * ((currentVolume - volume_) / dt_ - (volume_ - prevVolume_) / prevDt_));

                break;

            case 3:
                {
                    TFloat a = 0;
                    TFloat b = bulkModulus_ * ((isovolumetricRelaxationCavityVolume_ - currentVolume) / isovolumetricRelaxationCavityVolume_);
                    return std::max(a,b);
                    break;
                }

            case 4:
                return(ventricularPressure_ - (d_*(ventricularPressure_-ventricularResidualPressure_) * dt_) - (currentVolume -volume_) / c2_);

                break;

            default:
                throw std::runtime_error("void CBWindkessel3::CalcPressure(TFloat currentVolume): unsupported valve state");
        }


    else
        return(0);
}


TFloat CBWindkessel3::CalcPressureDerivative(TFloat vol)
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
                return( -b_  / initialVolume_ );

                break;

            case 0:
                return(0);

                break;

            case 1:
                return( -b_  / initialVolume_ );

                break;

            case 2:
                return(-(1+r1_/r2_) / c_ - r1_/dt_);

                break;

            case 3:
                return( -bulkModulus_  / isovolumetricRelaxationCavityVolume_ );

                break;

            case 4:
                return(-1 / c2_);

                break;

            default:
                throw std::runtime_error("void CBWindkessel3::CalcPressure(TFloat currentVolume): unsupported valve state");
        }


    else
        return(0);
}


CBWindkessel3::CBWindkessel3(CBCavity* cavity, ParameterMap* parameters) : CBCirculatorySystemModel(cavity, parameters)
{
    std::string index = std::to_string(Base::cavity_->GetIndex());

    Base::isActive_ = parameters->Get<bool>("Plugins.CirculatorySystem.Cavities.Cavity_" + index + ".Active", true);

    atrialPressure_ = parameters->Get<TFloat>("Plugins.CirculatorySystem.Cavities.Cavity_" + index + ".ClosingPressure", 2000);
    ventricularResidualPressure_ = parameters->Get<TFloat>("Plugins.CirculatorySystem.Cavities.Cavity_" + index + ".VentricularResidualPressure", 750);
    minSystolicPressure_         = parameters->Get<TFloat>("Plugins.CirculatorySystem.Cavities.Cavity_" + index + ".MinSystolicPressure", 10000.0);
    beatInterval_                = parameters->Get<TFloat>("Plugins.CirculatorySystem.BeatInterval", 1.);

    bulkModulus_ = parameters->Get<TFloat>("Plugins.CirculatorySystem.Cavities.Cavity_" + index + ".Windkessel3.BulkModulus", 1e7);
    c_           = parameters->Get<TFloat>("Plugins.CirculatorySystem.Cavities.Cavity_" + index + ".Windkessel3.C", 1e-8);
    c2_          = parameters->Get<TFloat>("Plugins.CirculatorySystem.Cavities.Cavity_" + index + ".Windkessel3.C2", 1e-7);
    r1_          = parameters->Get<TFloat>("Plugins.CirculatorySystem.Cavities.Cavity_" + index + ".Windkessel3.R1", 1.4e8);
    r2_          = parameters->Get<TFloat>("Plugins.CirculatorySystem.Cavities.Cavity_" + index + ".Windkessel3.R2", 1.4e8);
    d_           = parameters->Get<TFloat>("Plugins.CirculatorySystem.Cavities.Cavity_" + index + ".Windkessel3.D", 12);

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
        Base::file_ << "Time \t\t Volume \t\t  Ventricular_Pressure \t\t Aortic_Pressure \t\t dV/dt \t\t dPv/dt \t\t dPa/dt \t\t Volume-Pressure_Work \t\t Valve_State\n ";
        Base::file_.close();
    }
}
