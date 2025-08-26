//
//  CBRegistrationICS.cpp
//  CardioMechanics
//
//  Created by Thomas Fritz and Emily Reid on 06.09.11. 
//  Copyright (c) 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
//

// This is an extension of the non rigid iterative closest point algorithm with a point to surface approach


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

#include "CBRegistrationICS.h"
#include "CBSolver.h"


void CBRegistrationICS::Init()
{
    std::string targetNodesFilename = Base::GetParameters()->Get<std::string>("Plugins.NonRigidICS.TargetNodes", "");
	std::string targetSurfacesFilename = Base::GetParameters()->Get<std::string>("Plugins.NonRigidICS.TargetSurfaces", "");
	std::string targetOutputFilename = Base::GetParameters()->Get<std::string>("Plugins.NonRigidICS.TargetOutput", "Targets.vtu");
	atol_ = Base::GetParameters()->Get<TFloat>("Plugins.NonRigidICS.Atol", 0.0015);
    
    LoadTargetSurfaces(targetNodesFilename,targetSurfacesFilename);
    WriteSurfaceToFile(targetOutputFilename);
    
    std::string sourceFilename = Base::GetParameters()->Get<std::string>("Plugins.NonRigidICS.SourceSurfaces", "");
            
    LoadSourceNodesIndices(sourceFilename);
    
    alpha_ = Base::GetParameters()->Get<TFloat>("Plugins.NonRigidICS.Alpha", 0);
    beta_ = Base::GetParameters()->Get<TFloat>("Plugins.NonRigidICS.Beta", 1e-4);
    minDist_ = Base::GetParameters()->Get<TFloat>("Plugins.NonRigidICS.MinDist", 1e-4);
}

void CBRegistrationICS::WriteSurfaceToFile(const std::string& filename)
{
    vtkSmartPointer<vtkUnstructuredGrid> mesh = vtkSmartPointer<vtkUnstructuredGrid>::New();
    
    vtkSmartPointer<vtkPoints> nodes = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkIntArray> mat = vtkSmartPointer<vtkIntArray>::New();
    mat->SetNumberOfComponents(1);
    mat->SetName("ICS_Index");
   
    int n=0;
    for(std::map<TInt, std::vector<Triangle<TFloat> > >::iterator it = targetSurfaces_.begin(); it != targetSurfaces_.end(); it++)
    {
        for(std::vector<Triangle<TFloat> >::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++)
        {
            vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
            for(int i=0; i < 3; i++)
            {
                nodes->InsertNextPoint(it2->GetNode(i).X(),it2->GetNode(i).Y(),it2->GetNode(i).Z());
                mat->InsertNextValue(it->first);
                triangle->GetPointIds()->SetId(i, n);
                n++;
            }
             mesh->InsertNextCell(triangle->GetCellType(), triangle->GetPointIds());
        }
    }
    
    mesh->SetPoints(nodes);
    mesh->GetPointData()->AddArray(mat);

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(filename.c_str()); 
#if VTK_MAJOR_VERSION > 5
    writer->SetInputData(mesh);
#else
    writer->SetInput(mesh);
#endif
    writer->Write();
}

void CBRegistrationICS::LoadTargetSurfaces(const std::string& nodesFilename, const std::string& surfacesFilename)
{
    
    TFloat unit = Base::GetParameters()->Get<TFloat>("Plugins.NonRigidICS.Unit", 1);
    
    std::ifstream nodesFile(nodesFilename.c_str());
    std::ifstream surFile(surfacesFilename.c_str());

    if(!nodesFile.good())
        throw std::runtime_error("void CBNonRigidICP::Init(): File  " + nodesFilename + " does not exist");
    if(!surFile.good())
        throw std::runtime_error("void CBNonRigidICP::Init(): File  " + surfacesFilename + " does not exist");
        
        
    std::string str;
    
    int numNodes = 0;

    int numAttr = 0;
    
    int numElements = 0;
    int nodesPerElement = 0;
    
    getline(nodesFile, str);
    std::stringstream ss(str);
    
    ss >> numNodes >> nodesPerElement >> numAttr;
    
    getline(surFile, str);
    ss.str(str);
    ss.clear();
    
    ss >> numElements >> nodesPerElement >> numAttr;
    
    std::vector<Vector3<TFloat> > targetNodes;
    
    while(getline(nodesFile, str))
    {
        std::stringstream ss2(str);
        TInt index=-1;
        //            
        TFloat x,y,z;
        
        ss2 >> index >> x >> y >> z;
        
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
        targetNodes.push_back(Vector3<TFloat>(x,y,z));
    }

    
    nodesFile.close();
    
    
    while(getline(surFile, str))
    {
        std::stringstream ss2(str);
        
        TFloat index,n1,n2,n3;
        
        ss2 >> index >> n1 >> n2 >> n3;
        
        n1--;
        n2--;
        n3--;
        
        Triangle<TFloat> tri(targetNodes.at(n1),targetNodes.at(n2),targetNodes.at(n3));
        
        std::vector<TInt> attr;
        
        for(int i=0; i < numAttr; i++)
        {
            TInt t;
            ss2 >> t;
            attr.push_back(t);
        }
        
        std::map<TInt, std::vector<Triangle<TFloat> > >::iterator it = targetSurfaces_.find(attr.at(1));
        
        
        if( it == targetSurfaces_.end())
        {
            std::vector<Triangle<TFloat> > sur; 
            sur.push_back(tri);
            targetSurfaces_.insert(std::pair<TInt, std::vector<Triangle<TFloat> > >(attr.at(1),sur));
        }
        else
        {
            it->second.push_back(tri);
        }
    }
		surFile.close();
}

void CBRegistrationICS::Apply(PetscScalar time)
{
    vtkSmartPointer<vtkPoints> meshNodes = vtkSmartPointer<vtkPoints>::New();
		
	vtkSmartPointer<vtkUnstructuredGrid> mesh = vtkSmartPointer<vtkUnstructuredGrid>::New();

    
    if(time == 0)
        status_ = CBStatus::REPEAT;
    alpha_ = time * beta_;
    TFloat avgDist = 0;
	int cnt=0;
    TFloat a=0;
    TFloat b =0;
    
    
    for(std::map<TInt, std::vector<TInt> > ::iterator it = sourceNodesCoordsLocalIndices_.begin(); it != sourceNodesCoordsLocalIndices_.end(); it++)
    {
        
        
        
        std::map<PetscInt, std::vector<Triangle<PetscScalar> > >::iterator ts = targetSurfaces_.find(it->first); //adapt to triangle
        
        std::vector<Triangle<PetscScalar> >* cs = &(closestTriangle_.find(it->first)->second); //ok
		std::vector<bool>* cp = &(canBeProjected_.find(it->first)->second); //ok
        
        for(TInt i=0; i < it->second.size() / 3; i++)
        {
            PetscScalar sourceNodesCoords[3]; //ok
            Base::GetAdapter()->GetNodesCoords(3, &(it->second[3*i]) , sourceNodesCoords); //ok
            Vector3<PetscScalar> sourceNode(sourceNodesCoords[0],sourceNodesCoords[1],sourceNodesCoords[2]); 
            PetscScalar minDist = std::numeric_limits<PetscScalar>::max();
            PetscScalar dist = 0;
			Vector3<TFloat> tn;
            for(std::vector<Triangle<PetscScalar> >::iterator it2 = ts->second.begin(); it2 != ts->second.end(); it2++)
            {    
				Vector3<TFloat> t=it2->GetCentroid();
				
                dist = (sourceNode-t).Norm();
                for(int j=0; j < 3; j++)
                {
                    TFloat d2 = (sourceNode - it2->GetNode(j)).Norm();
                    if(d2 < dist)
                    {
                        dist = d2;
                        t = it2->GetNode(j);
                        
                    }
                }
                
				if( dist < minDist)
				{	
					Vector3<TFloat> nv = it2->GetNormalVector();
					Vector3<TFloat> ip ;
					ip = it2->CalcIntersectionPoint(sourceNode,nv);
					
					if(it2->GetDistanceTo(sourceNode + dist*nv) > dist)
						nv *= -1;
					
					if(it2->IsPointWithinTriangle(ip))
					{
                        dist = it2->GetDistanceTo(sourceNode);
						minDist = dist;
						cs->at(i) = *it2;
                        tn = ip;
                        cp->at(i) = true;
					}		
					else
					{
						minDist = dist;
                        tn = t;
						cs->at(i) = *it2;
                        cp->at(i) = false;
					}
                }
			}
            if(cp->at(i))
                a++;
            else
                b++;
			avgDist += minDist;
            GetAdapter()->GetSolver()->GetModel()->nodesScalarData.SetData("ICS_dist", it->second[3*i]/3, minDist);
            cnt++;
            
            meshNodes->InsertNextPoint(tn.X(),tn.Y(),tn.Z());
			meshNodes->InsertNextPoint(sourceNode(0),sourceNode(1),sourceNode(2));
            
                          


			vtkSmartPointer<vtkLine> line0 = vtkSmartPointer<vtkLine>::New();
			line0->GetPointIds()->SetId(0,2*cnt); //the second 0 is the index of the Origin in the vtkPoints
			line0->GetPointIds()->SetId(1,2*cnt+1);
			mesh->InsertNextCell(line0->GetCellType(),line0->GetPointIds());

  		}
	}

    avgDist /= cnt;
    cout << "Average Distance:" << avgDist << a/(a+b) << std::endl;
    
    
    mesh->SetPoints(meshNodes);
	vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	std::string filename;
	std::stringstream ss;
	ss << time*100;
	filename = "Closest/closestNeighbor." + ss.str() + ".vtu";
	writer->SetFileName(filename.c_str());
#if VTK_MAJOR_VERSION > 5
	writer->SetInputData(mesh);
#else
	writer->SetInput(mesh);
#endif
	writer->Write();
    
    if(avgDist < atol_)
        status_ = CBStatus::DACCORD;
}


void CBRegistrationICS::LoadSourceNodesIndices(const std::string& filename)
{
    
    if(filename != "")
    {
        std::ifstream sourceNodesIndicesFile(filename.c_str());
        
        if(!sourceNodesIndicesFile.good())
            throw std::runtime_error("void CBRegistrationICS::LoadSourceNodesIndices(): File  " + filename + " does not exist");
        
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
                
                if(nodes.insert(n).second == true)
                {
                    GetAdapter()->GetSolver()->GetModel()->nodesScalarData.SetData("ICS", n, a.at(1));

                    // nicht getestet !!!
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
                            it->second.push_back(n);
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
            std::vector<Triangle<TFloat> > vec;
            std::vector<bool> vec2;
            vec.resize(it->second.size(), Triangle<TFloat>());
            vec2.resize(it->second.size(), false);
            closestTriangle_.insert(std::pair<PetscInt, std::vector<Triangle<TFloat> > >(it->first,vec));
            canBeProjected_.insert(std::pair<PetscInt, std::vector<bool> >(it->first,vec2));
        }
        
        sourceNodesIndicesFile.close();
    }
    
}





void CBRegistrationICS::ApplyToNodalForces()
{
    TFloat avgForce=0;
    TInt cnt=0;
    for(std::map<TInt, std::vector<TInt> > ::iterator it = sourceNodesCoordsLocalIndices_.begin(); it != sourceNodesCoordsLocalIndices_.end(); it++)
    {
        std::vector<Triangle<PetscScalar> >* cs = &(closestTriangle_.find(it->first)->second); 
        std::vector<bool>* cp = &(canBeProjected_.find(it->first)->second);
        
        for(TInt i=0; i < it->second.size() / 3; i++)
        {
            PetscScalar sourceNodesCoords[3];
            Base::GetAdapter()->GetNodesCoords(3, &(it->second[3*i]) , sourceNodesCoords);
            Vector3<PetscScalar> sourceNode(sourceNodesCoords[0],sourceNodesCoords[1],sourceNodesCoords[2]);

            Triangle<TFloat> t = cs->at(i);
            
            Vector3<PetscScalar> regForce;
            
            if(cp->at(i))
            {
                Vector3<TFloat> nv = t.GetNormalVector();
                Vector3<TFloat> ip ;
                TFloat dist = t.GetDistanceTo(sourceNode);
                ip = t.CalcIntersectionPoint(sourceNode,nv);
                
                if(t.GetDistanceTo(sourceNode + dist*nv) > dist)
                    nv *= -1;

                regForce = -alpha_ * (ip - sourceNode);

            }
            else
            {
                Vector3<TFloat> bc=t.GetCentroid();
                regForce = -alpha_ * (bc - sourceNode);
            }
            avgForce += regForce.Norm();
            cnt++;
            PetscScalar regForceComponents[3] = {regForce(0), regForce(1), regForce(2)};
            Base::GetAdapter()->AddNodalForcesComponents(3,&(it->second[3*i]), regForceComponents);
        }
     
    }
    avgRegForce_ = avgForce / cnt;
}

void CBRegistrationICS::ApplyToNodalForcesJacobian()
{
    
    for(std::map<TInt, std::vector<TInt> > ::iterator it = sourceNodesCoordsLocalIndices_.begin(); it != sourceNodesCoordsLocalIndices_.end(); it++)
    {
        PetscScalar val[9] = {alpha_, 0, 0,0, alpha_, 0, 0, 0, alpha_};
        
        std::vector<Triangle<PetscScalar> >* cs = &(closestTriangle_.find(it->first)->second);
        std::vector<bool>* cp = &(canBeProjected_.find(it->first)->second);
        
        for(TInt i=0; i < it->second.size() / 3; i++)
        {
            TFloat sourceNodesCoords[3];
            Base::GetAdapter()->GetNodesCoords(3, &(it->second[3*i]) , sourceNodesCoords);
            Vector3<TFloat> sourceNode(sourceNodesCoords[0],sourceNodesCoords[1],sourceNodesCoords[2]);
            
            Triangle<TFloat> t = cs->at(i);
            
            Vector3<PetscScalar> regForce;
            
            if(cp->at(i))
            {
                Vector3<TFloat> nv = t.GetNormalVector();
                Vector3<TFloat> ip ;
                TFloat dist = t.GetDistanceTo(sourceNode);
                ip = t.CalcIntersectionPoint(sourceNode,nv);
                
                if(t.GetDistanceTo(sourceNode + dist*nv) > dist)
                    nv *= -1;
                
                regForce = -alpha_ * (ip - sourceNode);
            }
            else
            {
                Vector3<TFloat> bc=t.GetCentroid();
                regForce = -alpha_ * (bc - sourceNode);
            }

            
            TFloat e = Base::GetAdapter()->GetFiniteDifferencesEpsilon();
            for(int j=0;j<3;j++)
            {
                Vector3<TFloat> dS = sourceNode;
                dS(j) += e;
                Vector3<TFloat> df;
                if(cp->at(i))
                {
                    Vector3<TFloat> nv = t.GetNormalVector();
                    Vector3<TFloat> ip ;
                    TFloat dist = t.GetDistanceTo(dS);
                    ip = t.CalcIntersectionPoint(dS,nv);
                    
                    if(t.GetDistanceTo(dS + dist*nv) > dist)
                        nv *= -1;
                    
                    df = (-alpha_ * (ip - dS) - regForce) / e;
                }
                else
                {
                    Vector3<TFloat> bc=t.GetCentroid();
                    df = (-alpha_ * (bc - dS) - regForce) / e;
                }
                
                
                val[0 + j] = df(0);
                val[3 + j] = df(1);
                val[6 + j] = df(2);

            }

            Base::GetAdapter()->AddNodalForcesJacobianEntries(3, &(it->second[3*i]), 3, &(it->second[3*i]), val);
            
        }

        
    }
}

void CBRegistrationICS::FixSourceNodes()
{
    TInt numNodes = Base::GetAdapter()->GetSolver()->GetNumberOfNodes();
    fix_ = new bool[3*numNodes];
    for(int i=0; i < 3*numNodes; i++)
    {
        fix_[i] = false;
        GetAdapter()->GetSolver()->GetModel()->nodesScalarData.SetData("Fixed",i/3,0);
    }
    for(std::map<TInt, std::vector<TInt> >::iterator it = sourceNodesCoordsLocalIndices_.begin(); it != sourceNodesCoordsLocalIndices_.end(); it++)
        for(std::vector<TInt>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++)
        {
            fix_[*it2] = true;
            GetAdapter()->GetSolver()->GetModel()->nodesScalarData.SetData("Fixed",(*it2)/3,1);
        }  
    Base::GetAdapter()->LinkNodesComponentsBoundaryConditionsGlobal(fix_);
}



