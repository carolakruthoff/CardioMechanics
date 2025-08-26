/*
 *  CBModelLoaderTetgen.cpp
 *  CardioMechanics
 *
 *  Created by Nick van Osta on October 5th 2016.
 *  Copyright 2016 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#include "CBModelLoaderVTK.h"
#include "sstream"
#include "cstring"
#include "DCCtrl.h"

CBModel* CBModelLoaderVTK::Load(std::string Tag)
{
    Tag_ = Tag;
    std::string vtkFile     = Base::parameters_->Get<std::string>(Tag+".VTK.File", "");
    std::string type        = Base::parameters_->Get<std::string>(Tag+".VTK.Type", "");
    TFloat unit             = Base::parameters_->Get<TFloat>(Tag+".VTK.Unit", 1.0);
    bool determineNeighbors = Base::parameters_->Get<bool>(Tag+".VTK.DetermineNeighbors", true);
    
    LoadVelocity            = Base::parameters_->Get<bool>(Tag+".VTK.LoadVelocity",     true);
    LoadAcceleration        = Base::parameters_->Get<bool>(Tag+".VTK.LoadAcceleration", true);
    LoadSurfaceID           = Base::parameters_->Get<bool>(Tag+".VTK.LoadSurfaceID",    true);
    LoadActiveStress        = Base::parameters_->Get<bool>(Tag+".VTK.LoadActiveStress", true);
    LoadFixation            = Base::parameters_->Get<bool>(Tag+".VTK.Fixation",         true);
    LoadMaterial            = Base::parameters_->Get<bool>(Tag+".VTK.Material",         true);
    LoadFibers              = Base::parameters_->Get<bool>(Tag+".VTK.LoadFibers",       true);

    return(Load(type, vtkFile,unit,determineNeighbors));
}


CBModel* CBModelLoaderVTK::Load(std::string type, std::string vtkFile, TFloat unit, bool determineNeighbors)
{
    if(type=="vtu"){
        DCCtrl::print << "\n\nLoad VTU file...\n";
        vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
        reader->SetFileName(vtkFile.c_str());
        reader->Update();
        
        vtkSmartPointer<vtkUnstructuredGrid> mesh = reader->GetOutput();
        
        DCCtrl::print << "Load Nodes...\n";
        LoadNodes(mesh, unit);
        DCCtrl::print << "Load Elements...\n";
        LoadElements(mesh);
        
        if(determineNeighbors)
            model_->DetermineSolidElementNeighbors();
    }
    else{
        throw std::runtime_error("CBModelLoaderVTK::LoadData() type '" + type + "' is not known. For now, only 'vtu' is implemented");
        return nullptr;
    }
    
    DCCtrl::print << "VTK model loaded\n";
    return(Base::model_);
}


void CBModelLoaderVTK::LoadNodes(vtkSmartPointer<vtkUnstructuredGrid> Mesh, TFloat unit){
  auto numberPoints = Mesh->GetNumberOfPoints();
  std::vector<Vector3<TFloat>> nodes;
  nodes.reserve(numberPoints);
  std::vector<TInt> boundaryConditions;
  boundaryConditions.resize(numberPoints);
  
    for(int i=0; i<numberPoints; i++){
        TFloat *node = Mesh->GetPoint(i);
        Vector3<TFloat> nodeVec = Vector3<TFloat>(node);
        nodeVec*=unit; 
    //Base::model_->InsertNextNode(nodeVec);
    nodes.push_back(nodeVec);
    }
    
    //----------------------
    // Fixation
    //----------------------
    if(LoadFixation){
        DCCtrl::debug << "Load Fixation...\n";
        
        bool existFix;
        vtkDataArray *Fixation = GetArray(Mesh, "Fixation", "point", existFix);
        
        if(existFix)
        {
            for(int i=0; i<Fixation->GetSize(); i++)
      boundaryConditions[i] = Fixation->GetVariantValue(i).ToInt();
//    model_->SetNodeBoundaryCondition(i, Fixation->GetVariantValue(i).ToInt());
        }
    }

  model_->SetBoundaryConditions(std::move(boundaryConditions));
  model_->SetNodes(std::move(nodes));
  
    //----------------------
    // Velocity
    //----------------------
    if(LoadVelocity){
        DCCtrl::debug << "Load Velocity...\n";
        
        bool existVel;
        vtkDataArray *Velocity = GetArray(Mesh, "Velocity", "point", existVel, false);
        
        if(existVel){
            std::vector<Vector3<TFloat>> vel;
            for(int i=0; i<Velocity->GetSize(); i++){
                Vector3<TFloat> v = Vector3<TFloat>( Velocity->GetTuple3(i) );
                vel.push_back(v);
            }
            model_->SetVelocity(vel);
        }
    }

    //----------------------
    // Acceleration
    //----------------------
    if(LoadAcceleration){
        DCCtrl::debug << "Load Acceleration...\n";

        bool existAcc;
        vtkDataArray *Acceleration = GetArray(Mesh, "Acceleration", "point", existAcc, false);
        
        if(existAcc){
            std::vector<Vector3<TFloat>> acc;
            for(int i=0; i<Acceleration->GetSize(); i++){
                Vector3<TFloat> a = Vector3<TFloat>( Acceleration->GetTuple3(i) );
                acc.push_back(a);
            }
            model_->SetAcceleration(acc);
        }
    }

    //----------------------
    // General Point Data
    //----------------------
    DCCtrl::debug << "Load General Point Data...\n";
    
    int numberOfArrays = Mesh->GetPointData()->GetNumberOfArrays();
    for(int i = 0; i<numberOfArrays; i++){
        int numberOfComponents = Mesh->GetPointData()->GetArray(i)->GetElementComponentSize();
        
        if(numberOfComponents==1){
            for(int j=0; j<Mesh->GetPointData()->GetArray(i)->GetSize(); j++){
        model_->nodesScalarData.SetData(Mesh->GetPointData()->GetArray(i)->GetName(), j, *Mesh->GetPointData()->GetArray(i)->GetTuple(j));
            }
        }
        else if(numberOfComponents==3){
            for(int j=0; j<Mesh->GetPointData()->GetArray(i)->GetSize(); j++){
                Vector3<TFloat> tuple = Vector3<TFloat>(Mesh->GetPointData()->GetArray(i)->GetTuple3(j));
        model_->nodesVectorData.SetData(Mesh->GetPointData()->GetArray(i)->GetName(), j, tuple);
            }
        }
    }
}


void CBModelLoaderVTK::LoadElements(vtkSmartPointer<vtkUnstructuredGrid> Mesh){
    int numberCells = Mesh->GetNumberOfCells();
  std::vector<CBElement*> newElements;
  
  DCCtrl::debug << "Number of cells: " << numberCells << "\n";
    for(int i=0; i<numberCells; i++){
        vtkCell* cell = Mesh->GetCell(i);
        int type = cell->GetCellType();
        
        //Add elements here !!!
        switch(type) {
            case VTK_TRIANGLE: // Linear Triangle
            {
                CBElement* element = Base::model_->GetElementFactory()->New("T3");
                element->SetIndex(i);
                element->SetNodeIndex(0, cell->GetPointIds()->GetId(0));
                element->SetNodeIndex(1, cell->GetPointIds()->GetId(1));
                element->SetNodeIndex(2, cell->GetPointIds()->GetId(2));
                newElements.push_back(element);
                break;
            }
            case VTK_TETRA: // Linear Tetrahedron
            {
                CBElement* element = Base::model_->GetElementFactory()->New("T4");
                element->SetIndex(i);
                element->SetNodeIndex(0, cell->GetPointIds()->GetId(0));
                element->SetNodeIndex(1, cell->GetPointIds()->GetId(1));
                element->SetNodeIndex(2, cell->GetPointIds()->GetId(2));
                element->SetNodeIndex(3, cell->GetPointIds()->GetId(3));
                newElements.push_back(element);
                break;
            }
            case VTK_QUADRATIC_TRIANGLE: // Quadratic Triangle
            {
                CBElement* element = Base::model_->GetElementFactory()->New("T6");
                element->SetIndex(i);
                element->SetNodeIndex(0, cell->GetPointIds()->GetId(0));
                element->SetNodeIndex(1, cell->GetPointIds()->GetId(1));
                element->SetNodeIndex(2, cell->GetPointIds()->GetId(2));
                element->SetNodeIndex(3, cell->GetPointIds()->GetId(3));
                element->SetNodeIndex(4, cell->GetPointIds()->GetId(4));
                element->SetNodeIndex(5, cell->GetPointIds()->GetId(5));
                newElements.push_back(element);
                break;
            }
            case VTK_QUADRATIC_TETRA: // Quadratic Tetrahedron
            {
                CBElement* element = Base::model_->GetElementFactory()->New("T10");
                element->SetIndex(i);
                element->SetNodeIndex(0, cell->GetPointIds()->GetId(0));
                element->SetNodeIndex(1, cell->GetPointIds()->GetId(1));
                element->SetNodeIndex(2, cell->GetPointIds()->GetId(2));
                element->SetNodeIndex(3, cell->GetPointIds()->GetId(3));
                element->SetNodeIndex(4, cell->GetPointIds()->GetId(4));
                element->SetNodeIndex(5, cell->GetPointIds()->GetId(5));
                element->SetNodeIndex(6, cell->GetPointIds()->GetId(6));
                element->SetNodeIndex(7, cell->GetPointIds()->GetId(7));
                element->SetNodeIndex(8, cell->GetPointIds()->GetId(8));
                element->SetNodeIndex(9, cell->GetPointIds()->GetId(9));
                newElements.push_back(element);
                break;
            }
            default:
            {
                DCCtrl::print << "Error: element type " << type << " (" << vtkCellTypes::GetClassNameFromTypeId(type) << ") is not implemented in CBModelLoaderVTK::LoadElements() \n";
                throw std::runtime_error("Error: CBModelLoaderVTK::LoadElements()");
                break;
            }
        }
    }
    DCCtrl::debug << "ELEMENTS LOADED\n";
    
    //----------------------
    // Material
    //----------------------
    if(LoadMaterial){
        DCCtrl::debug << "Load Material...\n";
        
        bool existMat;
        vtkDataArray *Material = GetArray(Mesh, "Material", "cell", existMat, true);
        
        if(existMat)
        {
            for(int i=0; i<Material->GetSize(); i++){
      CBElement* element = model_->GetElements().at(i);
                element->SetMaterialIndex(Material->GetVariantValue(i).ToInt());
            }
        }
    }
    
    //----------------------
    // Surface Index
    //----------------------
    if(LoadSurfaceID){
        DCCtrl::debug << "Load Surface Index...\n";
        
        bool existSurEleIds, existSurIds;
        vtkDataArray *SurfaceElementIds = GetArray(Mesh, "SurfaceElementID", "cell", existSurEleIds, true);
        vtkDataArray *SurfaceIds = GetArray(Mesh, "SurfaceID", "cell", existSurIds, true);
        
        if(existSurEleIds && existSurIds)
        {
            for(int i=0; i<SurfaceElementIds->GetSize(); i++){
      CBElement* element = model_->GetElements().at(i);
                if(dynamic_cast<CBElementSurface*>(element)!=0){
                    CBElementSurface* s = dynamic_cast<CBElementSurface*>(element);
                    s->SetSurfaceElementIndex(SurfaceElementIds->GetVariantValue(i).ToInt());
                    s->SetSurfaceIndex(SurfaceIds->GetVariantValue(i).ToInt());
                }
            }
        }
    }
    
    //----------------------
    // Fibre Orientation
    //----------------------
    if(LoadFibers){
        DCCtrl::debug << "Load Fibre Orientation...\n";
        
        const int nQ = 5;
        bool ExistFiber[nQ];
        bool ExistSheet[nQ];
        bool ExistNormal[nQ];
        vtkDataArray *Fiber[nQ];
        vtkDataArray *Sheet[nQ];
        vtkDataArray *Normal[nQ];
        
        for(int i=0; i<nQ; i++){
            Fiber[i]  = GetArray(Mesh, ("Fiber"  + std::to_string(i)).c_str(), "cell", ExistFiber[i], true);
            Sheet[i]  = GetArray(Mesh, ("Sheet"  + std::to_string(i)).c_str(), "cell", ExistSheet[i], true);
            Normal[i] = GetArray(Mesh, ("Normal" + std::to_string(i)).c_str(), "cell", ExistNormal[i], true);
        }
        
    for(int i=0; i<model_->GetElements().size(); i++){
      CBElement* element = model_->GetElements().at(i);
            if(dynamic_cast<CBElementSolid*>(element) != 0){
                for(int j=0; j<element->GetNumberOfQuadraturePoints(); j++){
                    Vector3<TFloat> f;
                    Vector3<TFloat> s;
                    Vector3<TFloat> n;
                    
                    if(ExistFiber[j]){
                        double d[3];
                        Fiber[j]->GetTuple(i, d);
                        f = Vector3<TFloat>(d);
                    }
                    else
                        f = Vector3<TFloat>(1,0,0);
                    
                    if(ExistSheet[j])
                        s = Vector3<TFloat>(Sheet[j]->GetTuple3(i));
                    else
                        s =  Vector3<TFloat>(0,1,0);
                    
                    if(ExistNormal[j])
                        n = Vector3<TFloat>(Normal[j]->GetTuple3(i));
                    else
                        n =  Vector3<TFloat>(0,0,1);
                    
                    Matrix3<TFloat> b = Matrix3<TFloat>(f,s,n);
                    
                    element->SetBasisAtQuadraturePoint(j, b);
                }
            }
        }
    }
    
    //----------------------
    // ActiveStress
    //----------------------
    if(LoadActiveStress){
        DCCtrl::debug << "Load Active Stress...\n";
        
        bool existAS;
        vtkDataArray *ActiveStress = GetArray(Mesh, "ActiveStress", "cell", existAS, false);
        
        if(existAS){
            std::vector<TFloat> as;
            for(int i=0; i<ActiveStress->GetSize(); i++){
                TFloat a;
                ActiveStress->GetTuple(i, &a);
                as.push_back(a);
            }
            model_->SetActiveStress(as);
        }
    }
    
    //----------------------
    // General Cell Data
    //----------------------
    DCCtrl::debug << "Load General Cell Data...\n";
    
    int numberOfArrays = Mesh->GetCellData()->GetNumberOfArrays();
    for(int i = 0; i<numberOfArrays; i++){
        int numberOfComponents = Mesh->GetCellData()->GetArray(i)->GetElementComponentSize();
        
        if(numberOfComponents==1){
            for(int j=0; j<Mesh->GetCellData()->GetArray(i)->GetSize(); j++){
        model_->elementsScalarData.SetData(Mesh->GetCellData()->GetArray(i)->GetName(), j, *Mesh->GetCellData()->GetArray(i)->GetTuple(j));
            }
        }
        else if(numberOfComponents==3){
            for(int j=0; j<Mesh->GetCellData()->GetArray(i)->GetSize(); j++){
                Vector3<TFloat> tuple = Vector3<TFloat>(Mesh->GetCellData()->GetArray(i)->GetTuple3(j));
        model_->elementsVectorData.SetData(Mesh->GetCellData()->GetArray(i)->GetName(), j, tuple);
            }
        }
    }

  model_->SetElements(newElements);
}


vtkDataArray* CBModelLoaderVTK::GetArray(vtkSmartPointer<vtkUnstructuredGrid> Mesh, const char* Name, std::string type, bool &exists, bool required){
    vtkDataArray* A;
    exists = false;
    
    if(type=="cell"){
        if(Mesh->GetCellData()->HasArray(Name)){
            A = Mesh->GetCellData()->GetArray(Name);
            exists=true;
        }
        else if(required)
            DCCtrl::print << "CBModelLoaderVTK::LoadElements() The Cell array '" << Name << "' is not found in the source file. \n";
    }
    else if(type=="point"){
        if(Mesh->GetPointData()->HasArray(Name)){
            A = Mesh->GetPointData()->GetArray(Name);
            exists=true;
        }
        else if(required)
            DCCtrl::print << "CBModelLoaderVTK::LoadElements() The Point array '" << Name << "' is not found in the source file. \n";
    }
    else
        throw std::runtime_error("CBModelLoaderVTK::LoadElements() Type " + type + " is not known. Choose between 'cell' and 'point'.");
    
    if(required && !exists){
        vtkCellData *cd = Mesh->GetCellData();
        if(cd)
        {
            std::cout << " contains cell data with " << cd->GetNumberOfArrays() << " arrays." << std::endl;
            for (int i = 0; i < cd->GetNumberOfArrays(); i++)
            {
                std::cout << "\tArray " << i << " is named "
                << (cd->GetArrayName(i) ? cd->GetArrayName(i) : "NULL")
                << std::endl;
            }
        }
        throw std::runtime_error("CBModelLoaderVTK::LoadElements() The array is not found in the source file. ");
    }
    
    return A;
}

