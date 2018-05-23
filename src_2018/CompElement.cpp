//
//  GeoMesh.cpp
//  FemSC
//
//  Created by Philippe Devloo on 16/04/18.
//

#include "CompElement.h"
#include "CompElementTemplate.h"
#include "GeoElementTemplate.h"
#include "DataTypes.h"
#include "MathStatement.h"
#include "tpanic.h"

    // Default constructor of CompElement
CompElement::CompElement(){
    std::cout<<"bingo";
}
// Constructor of CompElement
CompElement::CompElement(int64_t ind, CompMesh *cmesh, GeoElement *geo){
    index = ind;
    compmesh = cmesh;
    geoel = geo;
    
}

// Copy constructor of CompElement
CompElement::CompElement(const CompElement &copy){
    index = copy.index;
    compmesh = copy.compmesh;
    geoel = copy.geoel;
    intrule = copy.intrule;
    mat = copy.mat;
    
}

// Operator of copy
CompElement &CompElement::operator=(const CompElement &copy){
   
    index = copy.index;
    compmesh = copy.compmesh;
    geoel = copy.geoel;
    intrule = copy.intrule;
    mat = copy.mat;
    return *this;

}

// Destructor of CompElement
CompElement::~CompElement(){
    
}

// Method for creating a copy of the element
CompElement *CompElement::Clone() const{
    
}

// Return the material object associated with the element
MathStatement *CompElement::GetStatement() const{
    return mat;
}

// Set the material object associated with the element
void CompElement::SetStatement(MathStatement *statement){
    mat=statement;
}

// Return integration rule established
IntRule *CompElement::GetIntRule() const{
    return intrule;
}

// Set integration rule established
void CompElement::SetIntRule(IntRule *intrul){
    intrule =intrul;
}


// Set element index
void CompElement::SetIndex(int64_t ind){
    index=ind;
}

// Return the geometric element associated
GeoElement *CompElement::GetGeoElement() const{
    return geoel;
}

// Set the geometric element associated
void CompElement::SetGeoElement(GeoElement *element){
    geoel = element;
}

// Return a pointer to the element computational mesh
CompMesh *CompElement::GetCompMesh() const{
    return compmesh;
    
}

// Set a pointer to the element computational mesh
void CompElement::SetCompMesh(CompMesh *mesh){
    compmesh = mesh;
}

// Initialize integration points data object
void CompElement::InitializeIntPointData(IntPointData &data) const{
    
  //revisar
    int dim = this->Dimension();
    int nshape = this->NShapeFunctions()+1; //quitar el mas 1
    const int nstate = this->GetStatement()->NState();
    //revisar.
    data.weight=0;
    data.x.resize(3);
    data.phi.resize(nshape,1);
    data.dphidx.Resize(dim,nshape);
    data.dphidksi.Resize(dim,nshape);
    data.axes.Resize(dim,3);
    data.gradx.Resize(dim, nshape);
    
    data.solution.resize(nstate);
    data.dsoldksi.Resize(dim,nstate);
    data.dsoldx.Resize(dim,nstate);
    
}

// Compute and fill integration points data object
void CompElement::ComputeRequiredData(IntPointData &data, VecDouble &intpoint) const{
    
    if (!geoel){
        DebugStop();
        std::cout << "CompEl ComputeRequiredData: no hay elemento asociado" << std::endl;
    }
  
    
    Matrix jac(1,1);
    data.ksi = intpoint;
    geoel->X(data.ksi, data.x);
    geoel->GradX(data.ksi, data.x, data.gradx);
    //this->ShapeFunctions(data.ksi, data.phi, data.dphidksi);
    
    //calcular axes y DetJac
    
    
    
    
    
    
    
}

// Compute the element stifness matrix and force vector
void CompElement::CalcStiff(Matrix &ek, Matrix &ef) const{
    
    MathStatement *material = mat;
    int nIntPoints = intrule->NPoints();
    IntPointData data;
    VecDouble intPoint;
    double weight;
    for(int i=0; i< nIntPoints;i++){
        //inicializa los datos a utilizar
        this->InitializeIntPointData(data);
        //devuelve el punto de integracion y el peso
        intrule->Point(i, intPoint , weight);
        this->ComputeRequiredData(data, intPoint);
        //devuelve el punto de integracion y el peso
        weight *=fabs(data.detjac);
        material->Contribute(data, weight, ek, ef);
        
    }

}

