//
//  CBApplyPressureFromFile.h
//  CardioMechanics
//
//  Created by ek717 on 24.04.20.
//
//

#ifndef CB_APPLY_PRESSURE_FROM_FILE
#define CB_APPLY_PRESSURE_FROM_FILE

#include <map>
#include <memory>

#include "CBSolverPlugin.h"
#include "CBDataFromPressureFile.h"
#include "CBCirculationCavity.h"
#include "CBSolverPluginCavities.h"

class CBApplyPressureFromFile : public CBSolverPluginCavities
{
public:
    CBApplyPressureFromFile(){}
    virtual ~CBApplyPressureFromFile(){}
    
    virtual void Init() override;
    virtual void Apply(TFloat time) override;
    virtual void ApplyToNodalForces() override;
    virtual void ApplyToNodalForcesJacobian() override;
    virtual void AnalyzeResults() override;
    virtual void WriteHeaderToFile();
    virtual void WriteToFile(TFloat time) override;
    bool WantsToAnalyzeResults() { return true; }
    CBStatus GetStatus() override {return status_;}
    virtual std::string GetName() override { return("CBApplyPressureFromFile"); };
    
private:
    std::ofstream fileIn_;
    std::string filenameIn_;
    std::ofstream fileOut_;
    std::string filenameOut_;
    bool headerWritten_ = false;
    
    struct IntervalStruct
    {
        std::string pressureFileName;
        std::shared_ptr<CBDataFromPressureFile> pressureFile;
        bool relaxElementsAtStart;
        std::vector<TInt> materialsToRelax;
        bool invert;
        TFloat startTime;
        TFloat stopTime;
        TFloat offset;
        TFloat amplitude;
    };
    
    struct ValueStruct
    {
        TFloat pressure;
        TFloat volume;
    };
    
    std::map<TInt, std::vector<IntervalStruct>> groups_;
    std::map<TInt, std::vector<IntervalStruct>>::iterator groupsIt_;
    std::map<TInt, ValueStruct> valuesStructs_;
    std::map<TInt, ValueStruct>::iterator valuesIt_;
    
    std::map<TInt, TFloat> materialsRelaxed_;
    std::map<TInt, TFloat>::iterator materialsRelaxedIt_;
    
    typedef CBSolverPlugin Base;
};

#endif
