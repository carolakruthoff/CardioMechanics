//
//  CBPointsCtrl.h
//  CardioMechanics_unstable
//
//  Created by Thomas Fritz on 09.03.13.
//
//

#ifndef CB_POINTS_CTRL_H
#define CB_POINTS_CTRL_H

#include "CBSolverPlugin.h"
#include "CBDataFromFile.h"

class CBPointsCtrl : public CBSolverPlugin
{
public:
    ~CBPointsCtrl(){if(pointsCoords_)delete pointsCoords_;}
    virtual void Apply(TFloat time);
    
    virtual std::string GetName(){return(std::string("PointsControl"));}
    virtual void Init();

protected:
    CBDataFromFile*   pointsCoords_;
    std::vector<TInt> controllNodesLocalIndices_;
private:
};
#endif