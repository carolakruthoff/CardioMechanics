#pragma once

#include <string>
#include <vector>

class CBElement;
class CBSolver;
class CBConstitutiveModel;
class CBCirculatorySystemModel;
class ParameterMap;
class CBSolverPlugin;
class CBPreprocessing;
class CBCavity;


CBElement* CBExtensionExperimentalNewElement(std::string elementType);
CBSolver* CBExtensionExperimentalLoadSolver(std::string solverType);
CBConstitutiveModel* CBExtensionExperimentalLoadConstitutiveModel(std::string solverType);
CBCirculatorySystemModel* CBExtensionExperimentalLoadCirculatorySystemModel(CBCavity* cavity, ParameterMap* parameters);
void CBExtensionExperimentalLoadPreprocessingPlugins(std::vector<CBPreprocessing* >& plugins, ParameterMap* parameters);
void CBExtensionExperimentalLoadSolverPlugins(std::vector<CBSolverPlugin* >& plugins, ParameterMap* parameters);
