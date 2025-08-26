/*
 *  CBModelLoaderTetgen.h
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 13.04.10.
 *  Copyright 2010 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */
#ifndef CB_MODEL_LOADER_VTK
#define CB_MODEL_LOADER_VTK

#include <iostream>
#include <fstream>
#include <string>

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
#include <vtkQuadraticTriangle.h>

#include <vtkUnstructuredGrid.h>
#include <vtkIntArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkProperty.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <vtkGenericDataObjectReader.h>
#include <vtkXMLUnstructuredGridReader.h>

#include "CBModelLoader.h"


class CBModelLoaderVTK : public CBModelLoader
{
public:
    CBModelLoaderVTK(ParameterMap* parameters) : Base(parameters){};
    virtual ~CBModelLoaderVTK(){}
    CBModel* Load(std::string Tag = "Mesh");
    CBModel* Load(std::string type,std::string vtkFile, TFloat unit, bool determineNeighbors=true);
protected:
    vtkDataObject* LoadData(std::string type, std::string vtkFile, TFloat unit);
    void LoadNodes(vtkSmartPointer<vtkUnstructuredGrid> Mesh, TFloat unit);
    void LoadElements(vtkSmartPointer<vtkUnstructuredGrid> Mesh);
    vtkDataArray* GetArray(vtkSmartPointer<vtkUnstructuredGrid> Mesh, const char* Name, std::string type, bool &exists, bool required = true);
    
    
    vtkDataObject* Data;
    
    bool IsDataLoaded = false;
private:
    typedef CBModelLoader   Base;
    std::string Tag_;
    
    bool LoadVelocity=true;
    bool LoadAcceleration=true;
    bool LoadSurfaceID=true;
    bool LoadActiveStress=true;
    bool LoadFixation=true;
    bool LoadMaterial=true;
    bool LoadFibers = true;
};

#endif
