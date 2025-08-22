#pragma once
#include <string>
#include <vector>

class ParameterMap;
class CBElement;
class CBSolver;
class CBConstitutiveModel;
class CBCirculatorySystemModel;
class CBCavity;
class CBSolverPlugin;
class CBPreprocessing;

CBElement* CBExtensionLoadElement(std::string elementType);
CBSolver* CBExtensionLoadSolver(std::string solverType);
CBConstitutiveModel* CBExtensionLoadConstitutiveModel(std::string solverType);
CBCirculatorySystemModel* CBExtensionLoadCirculatorySystemModel(CBCavity* cavity, ParameterMap* parameters);
void CBExtensionLoadSolverPlugins(std::vector<CBSolverPlugin* >& plugins, ParameterMap* parameters);
void CBExtensionLoadPreprocessingPlugins(std::vector<CBPreprocessing* >& plugins, ParameterMap* parameters);
