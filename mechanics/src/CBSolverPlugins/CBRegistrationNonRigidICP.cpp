//
//  CBRegistrationNonRigidICP.cpp
//  CardioMechanics
//
//  Created by Thomas Fritz on 24.11.11.
//  Copyright (c) 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
//

#include <fstream>
#include <string>
#include <sstream>

#include <iostream>
#include <limits>

#include "Matrix3.h"

#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkLine.h>
#include <vtkTriangle.h>
#include <vtkPolyData.h>
#include <vtkDataSetMapper.h>
#include <vtkPlane.h>
#include <vtkTetra.h>
#include <vtkQuadraticTetra.h>
#include <vtkUnstructuredGrid.h>
#include <vtkIntArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkProperty.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include "CBRegistrationNonRigidICP.h"
#include "CBSolver.h"


void CBRegistrationNonRigidICP::Init()
{
    std::string targetFile = Base::GetParameters()->Get<std::string>("Plugins.NonRigidICP.TargetNodes", "");
    LoadTargetNodes(targetFile);
    std::string sourceFile = Base::GetParameters()->Get<std::string>("Plugins.NonRigidICP.SourceNodesIndices", "");
    LoadSourceNodesIndices(sourceFile);
    alpha_ = Base::GetParameters()->Get<TFloat>("Plugins.NonRigidICP.Alpha", 0);
    beta_ = Base::GetParameters()->Get<TFloat>("Plugins.NonRigidICP.Beta", 1e-4);
    minDist_ = Base::GetParameters()->Get<TFloat>("Plugins.NonRigidICP.MinDist", 1e-4);
}


void CBRegistrationNonRigidICP::LoadTargetNodes(const std::string& filename)
{

    TFloat unit = Base::GetParameters()->Get<TFloat>("Plugins.NonRigidICP.Unit", 1);
    int n=0;
    if(filename != "")
    {
        std::ifstream targetNodesFile(filename.c_str());
        
        if(!targetNodesFile.good())
            throw std::runtime_error("void CBNonRigidICP::Init(): File  " + filename + " does not exist");
        
        std::string str;
        
        
        int numNodes = 0;
        int nodesPerElement = 0;
        int numAttr = 0;
        
        getline(targetNodesFile, str);
        std::stringstream ss(str);
        
        ss >> numNodes >> nodesPerElement >> numAttr;

        vtkSmartPointer<vtkPoints> meshNodes = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkIntArray> mat = vtkSmartPointer<vtkIntArray>::New();
        mat->SetNumberOfComponents(1);
        mat->SetName("Index");
        
        while(getline(targetNodesFile, str))
        {
            std::stringstream ss2(str);
            TInt index=-1;
            
            TFloat x,y,z;
            
            ss2 >> index >> x >> y >> z;

            if(index ==-1)
                continue;
            else 
                n++;
            
            x *= unit;
            y *= unit;
            z *= unit; 
            
            std::vector<TInt> a;
            
            for(int i=0; i < numAttr; i++)
            {
                TInt t;
                ss2 >> t;
                a.push_back(t);
            }

            std::map<TInt, std::vector<Vector3<TFloat> > >::iterator it = targetNodes_.find(a.at(1));
            if( it == targetNodes_.end())
            {
                std::vector<Vector3<TFloat> > nodes;
                nodes.push_back(Vector3<TFloat>(x,y,z));
                targetNodes_.insert(std::pair<TInt, std::vector<Vector3<TFloat> > >(a.at(1),nodes));
            }
            else
            {
                it->second.push_back(Vector3<TFloat>(x,y,z));
            }
            mat->InsertNextValue(a.at(1));
            meshNodes->InsertNextPoint(x,y,z);            
        }

        if(numNodes != n)
            throw  std::runtime_error("void CBRegistrationNonRigidICP::LoadTargetNodes(): File " + filename + " is corrupt. Number of nodes does not fit to number of lines ");
        
        targetNodesFile.close();

        
        vtkSmartPointer<vtkUnstructuredGrid> mesh = vtkSmartPointer<vtkUnstructuredGrid>::New();
        mesh->SetPoints(meshNodes);
        mesh->GetPointData()->AddArray(mat);
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        writer->SetFileName("TargetNodes.vtu"); 
#if VTK_MAJOR_VERSION > 5
        writer->SetInputData(mesh);
#else
        writer->SetInput(mesh);
#endif
        writer->Write();

    }
}






void CBRegistrationNonRigidICP::LoadSourceNodesIndices(const std::string& filename)
{
    
    if(filename != "")
    {
        std::ifstream sourceNodesIndicesFile(filename.c_str());
        
        if(!sourceNodesIndicesFile.good())
            throw std::runtime_error("void CBRegistrationNonRigidICP::LoadSourceNodesIndices(): File  " + filename + " does not exist");
        
        std::string str;
        std::map<PetscInt, std::vector<PetscInt> > indices;
        
        int numNodes = 0;
        int nodesPerElement = 0;
        int numAttr = 0;
        
        getline(sourceNodesIndicesFile, str);
        std::stringstream ss(str);
        
        ss >> numNodes >> nodesPerElement >> numAttr;
        
        std::set<int> nodes;    
        
        while(getline(sourceNodesIndicesFile, str))
        {
            std::stringstream ss2(str);
            PetscInt n;
            ss2 >> n;
            std::vector<int> v;
            std::vector<int> a;
            
            for(int i=0; i < nodesPerElement; i++)
            {
                PetscInt t;
                ss2 >> t;
                v.push_back(t-1);
                
            }
            
            for(int i=0; i < numAttr; i++)
            {
                PetscInt t;
                ss2 >> t;
                a.push_back(t);
            }
            
            for(int i=0; i < v.size(); i++)
            {
                
                n = v.at(i);
                
                n = GetAdapter()->GetSolver()->GetModel()->GetForwardMapping(n);
                
                if(nodes.insert(n).second == true && a.at(1) == 3)
                {
                    GetAdapter()->GetSolver()->GetModel()->nodesScalarData.SetData("ICP", n, a.at(1));
                    //                     if(n >= Base::GetAdapter()->GetLocalNodesFrom() && n <= Base::GetAdapter()->GetLocalNodesTo() )
//#warning nicht getestet !!!
                    if(n >= Base::GetAdapter()->GetSolver()->GetLocalNodesFrom() && n <= Base::GetAdapter()->GetSolver()->GetLocalNodesTo() )
                    {
                        std::map<PetscInt, std::vector<PetscInt> >::iterator it = indices.find(a.at(1));
                        if( it == indices.end())
                        {
                            std::vector<PetscInt> vec;
                            vec.push_back(n);
                            indices.insert(std::pair<PetscInt, std::vector<PetscInt> >(a.at(1),vec));
                        }
                        else
                        {
                            it->second.push_back(n);
                        }
                    }
                }
            }
            
        }
        
        
        
        for(std::map<PetscInt, std::vector<PetscInt> >::iterator it = indices.begin(); it != indices.end(); it++)
        {

            std::vector<PetscInt> nodesCoordsIndices;
            nodesCoordsIndices.resize(3 * it->second.size(),0);
            int i=0;
            for(std::vector<PetscInt>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++)
            {
                //                 PetscInt localIndex = *it2 - Base::GetAdapter()->GetLocalNodesFrom();
//#warning nicht getestet !!!!
                
                PetscInt localIndex = *it2 - Base::GetAdapter()->GetSolver()->GetLocalNodesFrom();

                nodesCoordsIndices[3*i    ] = 3* localIndex;
                nodesCoordsIndices[3*i + 1] = 3* localIndex + 1;
                nodesCoordsIndices[3*i + 2] = 3* localIndex + 2;
                i++;
            }
            
            sourceNodesCoordsLocalIndices_.insert(std::pair<TInt, std::vector<TInt> >(it->first, nodesCoordsIndices));
            std::vector<Vector3<PetscScalar> > vec;
            vec.resize(it->second.size(), Vector3<PetscScalar>(0,0,0));
            closestNeighbors_.insert(std::pair<PetscInt, std::vector<Vector3<PetscScalar> > >(it->first,vec));
        }
        
        sourceNodesIndicesFile.close();
    }
}

void CBRegistrationNonRigidICP::Apply(PetscScalar time)
{
    alpha_ = time * beta_;
    
    
    for(std::map<TInt, std::vector<TInt> > ::iterator it = sourceNodesCoordsLocalIndices_.begin(); it != sourceNodesCoordsLocalIndices_.end(); it++)
    {

        std::map<PetscInt, std::vector<Vector3<PetscScalar> > >::iterator it2 = targetNodes_.find(it->first);
        
        std::vector<Vector3<PetscScalar> >* cn = &(closestNeighbors_.find(it->first)->second);
        for(TInt i=0; i < it->second.size() / 3; i++)
        {
            PetscScalar sourceNodesCoords[3];
            Base::GetAdapter()->GetNodesCoords(3, &(it->second[3*i]) , sourceNodesCoords);
            Vector3<PetscScalar> sourceNode(sourceNodesCoords[0],sourceNodesCoords[1],sourceNodesCoords[2]);
            PetscScalar minDistance = std::numeric_limits<PetscScalar>::max();
            PetscScalar dist = 0;
            
            Vector3<PetscScalar> closestNeighbor;
            
            for(std::vector<Vector3<PetscScalar> >::iterator it3 = it2->second.begin(); it3 != it2->second.end(); it3++)
            {                
                dist = (*it3 - sourceNode).Norm();
                if(  dist < minDistance)
                {
                    closestNeighbor = *it3;
                    minDistance = dist;
                }
            }
            cn->at(i) = closestNeighbor;
        }
    }
    
}

void CBRegistrationNonRigidICP::ApplyToNodalForces()
{
    
    for(std::map<TInt, std::vector<TInt> > ::iterator it = sourceNodesCoordsLocalIndices_.begin(); it != sourceNodesCoordsLocalIndices_.end(); it++)
    {
  //  #warning Sehr bedenklich
        std::vector<Vector3<PetscScalar> >* cn = &(closestNeighbors_.find(it->first)->second);
        for(TInt i=0; i < it->second.size() / 3; i++)
        {
            PetscScalar sourceNodesCoords[3];
            Base::GetAdapter()->GetNodesCoords(3, &(it->second[3*i]) , sourceNodesCoords);            Vector3<PetscScalar> sourceNode(sourceNodesCoords[0],sourceNodesCoords[1],sourceNodesCoords[2]);
            Vector3<PetscScalar> regForce = -alpha_ * (cn->at(i) - sourceNode);
            PetscScalar regForceComponents[3] = {regForce(0), regForce(1), regForce(2)};
            Base::GetAdapter()->AddNodalForcesComponents(3,&(it->second[3*i]), regForceComponents);
        }
    }
}

void CBRegistrationNonRigidICP::ApplyToNodalForcesJacobian()
{
    
     for(std::map<TInt, std::vector<TInt> > ::iterator it = sourceNodesCoordsLocalIndices_.begin(); it != sourceNodesCoordsLocalIndices_.end(); it++)
     {
        for(PetscInt i=0; i < it->second.size() / 3; i++)
        {
            PetscScalar val[9] = {alpha_, 0, 0,0, alpha_, 0, 0, 0, alpha_};
            Base::GetAdapter()->AddNodalForcesJacobianEntries(3, &(it->second[3*i]), 3, &(it->second[3*i]), val);
        }
        
    }
}






