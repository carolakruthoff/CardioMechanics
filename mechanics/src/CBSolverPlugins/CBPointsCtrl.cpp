//
//  CBPointsCtrl.h.cpp
//  CardioMechanics_unstable
//
//  Created by Thomas Fritz on 09.03.13.
//
//

#include <vector>

#include "CBPointsCtrl.h"
#include "CBSolver.h"

void CBPointsCtrl::Init()
{
    std::string filename = parameters_->Get<std::string>("Plugins.PointsCtrl.Coordinates.File");

    pointsCoords_ = new CBDataFromFile();

    pointsCoords_->Init(filename);
    filename = parameters_->Get<std::string>("Plugins.PointsCtrl.Points.File","");
    bool* bc = GetAdapter()->GetSolver()->GetNodesComponentsBoundaryConditionsGlobal();
    if(filename != "")
    {
        std::ifstream controllNodesFile(filename.c_str());

        std::string str;

        if(!controllNodesFile.good())
            throw std::runtime_error("void CBPointsCtrl.h::Init(): File " + filename + " does not exist");

        while(getline(controllNodesFile, str))
        {
            std::stringstream ss(str);

            int n;
            ss >> n;
            n--;
            bc[3*n] = true;
            bc[3*n+1] = true;
            bc[3*n+2] = true;
            n = GetAdapter()->GetSolver()->GetModel()->GetForwardMapping(n);

            if((n >= GetAdapter()->GetSolver()->GetLocalNodesFrom()) && (n <= GetAdapter()->GetSolver()->GetLocalNodesTo()))
                controllNodesLocalIndices_.push_back(n - GetAdapter()->GetSolver()->GetLocalNodesFrom());
        }
    }
    else
    {
        std::vector<TInt> fromTo = parameters_->GetArray<TInt>("Plugins.PointsCtrl.Points.FromTo");

        if(fromTo.size() != 2)
            throw std::runtime_error("void BPointsCtrl.h::Init(): Plugins.PointsCtrl.Points.FromTo may only contain two values");

        for(int i = fromTo[0]-1; i <= fromTo[1]-1; i++)
        {
            bc[3*i] = true;
            bc[3*i+1] = true;
            bc[3*i+2] = true;
            
            TInt n = GetAdapter()->GetSolver()->GetModel()->GetForwardMapping(i);
            
            if((n >= GetAdapter()->GetSolver()->GetLocalNodesFrom()) && (n <= GetAdapter()->GetSolver()->GetLocalNodesTo()))
                controllNodesLocalIndices_.push_back(n - GetAdapter()->GetSolver()->GetLocalNodesFrom());
        }
    }

    GetAdapter()->LinkNodesComponentsBoundaryConditionsGlobal(bc);
}


void CBPointsCtrl::Apply(TFloat time)
{
    for(int i=0; i < controllNodesLocalIndices_.size(); i++)
    {
        TFloat newCoords[3] = {pointsCoords_->Get(time, 3*i), pointsCoords_->Get(time, 3*i+1), pointsCoords_->Get(time, 3*i+2)};
        TInt newCoordsIndices[3] = { 3*controllNodesLocalIndices_[i],3*controllNodesLocalIndices_[i]+1,3*controllNodesLocalIndices_[i]+2};
        GetAdapter()->SetNodesCoords(3,newCoordsIndices,newCoords);
    }
}