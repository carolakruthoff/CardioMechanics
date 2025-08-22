#include "CBExtensionExperimental.h"

#include "CBSolver.h"
#include "CBElement.h"
#include "CBCirculatorySystemModel.h"
#include "CBConstitutiveModel.h"

// --------------- Experimental Elements --------------
#include "CBElementSolidT10RI.h"
#include "CBElementSolidT10TH.h"
#include "CBElementSolidT10T4.h"
#include "CBElementSolidT10RIT4.h"
// -----------------------------------------------------

// ---------------- Experimental Solver ---------------
#include "CBSolverActiveStressEstimator.h"
#include "CBSolverActiveStressTensorEstimator.h"
#include "CBSolverActiveStressEstimatorNewmarkBeta.h"
// -----------------------------------------------------

// ------------ Experimental Solver Plugins -----------
#include "CBDetermineNodalForces.h"
#include "CBRegistration.h"
#include "CBRegistrationICS.h"
#include "CBRegistrationNonRigidICP.h"
// -----------------------------------------------------

// --------- Experimental Preprocessing Plugins --------
// #include "..."
// -----------------------------------------------------

// --------- Experimental Constitutive Model -----------
#include "CBConstitutiveModelCosta.h"
// -----------------------------------------------------

// ---------- Experimental Circulatory Model -----------
#include "CBPressureVolumeData.h"
#include "CBCirculatorySystemModelFromFile.h"
// -----------------------------------------------------

CBElement* CBExtensionExperimentalNewElement(std::string elementType)
{
    CBElement* element = 0;
    if(elementType == std::string("T10RI"))
        element = new CBElementSolidT10RI;
    else if(elementType == std::string("T10TH") || elementType == std::string("T10TH2") /*Backward Compability*/)
        element = new CBElementSolidT10TH;
    else if(elementType == std::string("T10T4"))
        element = new CBElementSolidT10T4;
    else if(elementType == std::string("T10RIT4"))
        element = new CBElementSolidT10RIT4;
    return element;
}

CBSolver* CBExtensionExperimentalLoadSolver(std::string solverType)
{
    if(solverType == "ActiveStressEstimator")
        return new CBSolverActiveStressEstimator();
    if(solverType == "ActiveStressTensorEstimator")
        return new CBSolverActiveStressTensorEstimator();
    else if(solverType == "ActiveStressEstimatorNewmarkBeta")
        return new CBSolverActiveStressEstimatorNewmarkBeta();
    else
        return 0;
}

void CBExtensionExperimentalLoadSolverPlugins(std::vector<CBSolverPlugin* >& plugins, ParameterMap* parameters)
{
    if(parameters->Get<bool>("Solver.Plugins.NonRigidICP", false))
    {
        CBSolverPlugin* solverPlugin = new CBRegistrationNonRigidICP();
        solverPlugin->SetParameters(parameters);
        plugins.push_back(solverPlugin);
    }
 
    if(parameters->Get<bool>("Solver.Plugins.DetermineNodalForces", false))
    {
        CBSolverPlugin* solverPlugin = new CBDetermineNodalForces();
        solverPlugin->SetParameters(parameters);
        plugins.push_back(solverPlugin);
    }
    
}

void CBExtensionExperimentalLoadPreprocessingPlugins(std::vector<CBPreprocessing* >& plugins, ParameterMap* parameters)
{
    
}

CBConstitutiveModel* CBExtensionExperimentalLoadConstitutiveModel(std::string modelType)
{
    if(modelType == "Costa")
        return new CBConstitutiveModelCosta();
    else
        return 0;
}

CBCirculatorySystemModel* CBExtensionExperimentalLoadCirculatorySystemModel(CBCavity* cavity, ParameterMap* parameters)
{
    std::string type = parameters->Get<std::string>("Plugins.CirculatorySystem.Cavities.Cavity_" + std::to_string(cavity->GetIndex()) + ".Type","");
    CBCirculatorySystemModel* model = 0;
    if(type == "PVDiagramm")
        model = new CBPressureVolumeData(cavity, parameters);
    if(type == "FromFile")
        model = new CBCirculatorySystemModelFromFile(cavity, parameters);
    return model;
}