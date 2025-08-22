#include "CBExtension.h"
#include <string>
#include <stdexcept>

#ifdef EXPERIMENTAL
#include "CBExtensionExperimental.h"
#endif

#ifdef CONTRIBUTIONS
#include "CBExtensionContribution.h"
#endif

CBElement* CBExtensionLoadElement(std::string elementType)
{
    CBElement* e1 = 0;
    CBElement* e2 = 0;
#ifdef EXPERIMENTAL
    e1 = CBExtensionExperimentalNewElement(elementType);
#endif
    
#ifdef CONTRIBUTIONS
    e2 = CBExtensionContributionNewElement(elementType);
#endif
    if(e1 && e2)
    {
        throw std::runtime_error("CBElement* CBExtensionLoadElement(std::string elementType): Element type is defined in EXPERIMENTAL and CONTRIBUTIONS");
    }
    else
    {
        if(e1)
            return e1;
        else if(e2)
            return e2;
        else
            return 0;
    }
}

CBSolver* CBExtensionLoadSolver(std::string solverType)
{
    CBSolver* s1 = 0;
    CBSolver* s2 = 0;
#ifdef EXPERIMENTAL
    s1 = CBExtensionExperimentalLoadSolver(solverType);
#endif
    
#ifdef CONTRIBUTIONS
    s2 = CBExtensionContributionLoadSolver(solverType);
#endif

    if(s1 && s2)
    {
        throw std::runtime_error("CBSolver* CBExtensionLoadSolver(std::string solverType): Solver type is defined in EXPERIMENTAL and in CONTRIBUTIONS");
    }
    else
    {
        if(s1)
            return s1;
        else if(s2)
            return s2;
        else
            return 0;
    }
}

CBConstitutiveModel* CBExtensionLoadConstitutiveModel(std::string modelType)
{
    CBConstitutiveModel* m1 = 0;
    CBConstitutiveModel* m2 = 0;
    
#ifdef EXPERIMENTAL
    m1 = CBExtensionExperimentalLoadConstitutiveModel(modelType);
#endif
    
#ifdef CONTRIBUTIONS
    m2 = CBExtensionContributionLoadConstitutiveModel(modelType);
#endif
    if(m1 && m2)
    {
        throw std::runtime_error("CBSolver* CBExtensionLoadSolver(std::string solverType): Solver type is defined in EXPERIMENTAL and in CONTRIBUTIONS");
    }
    else
    {
        if(m1)
            return m1;
        else if(m2)
            return m2;
        else
            return 0;
    }
}


CBCirculatorySystemModel* CBExtensionLoadCirculatorySystemModel(CBCavity* cavity, ParameterMap* parameters)
{
    
    CBCirculatorySystemModel* m1 = 0;
    CBCirculatorySystemModel* m2 = 0;
    
#ifdef EXPERIMENTAL
    m1 = CBExtensionExperimentalLoadCirculatorySystemModel(cavity, parameters);
#endif
    
#ifdef CONTRIBUTIONS
    m2 = CBExtensionContributionLoadCirculatorySystemModel(cavity, parameters);
#endif
    if(m1 && m2)
    {
        throw std::runtime_error("CBSolver* CBExtensionLoadSolver(std::string solverType): Solver type is defined in EXPERIMENTAL and in CONTRIBUTIONSS");
    }
    else
    {
        if(m1)
            return m1;
        else if(m2)
            return m2;
        else
            return 0;
    }
}

void CBExtensionLoadPreprocessingPlugins(std::vector<CBPreprocessing* >& plugins, ParameterMap* parameters)
{
    
#ifdef EXPERIMENTAL
    CBExtensionExperimentalLoadPreprocessingPlugins(plugins,parameters);
#endif
    
#ifdef CONTRIBUTIONS
    CBExtensionContributionLoadPreprocessingPlugins(plugins,parameters);
#endif
}

void CBExtensionLoadSolverPlugins(std::vector<CBSolverPlugin* >& plugins, ParameterMap* parameters)
{
    
#ifdef EXPERIMENTAL
    CBExtensionExperimentalLoadSolverPlugins(plugins,parameters);
#endif
    
#ifdef CONTRIBUTIONS
    CBExtensionContributionLoadSolverPlugins(plugins,parameters);
#endif
}
