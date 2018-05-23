//
//  GeoMesh.cpp
//  FemSC
//
//  Created by Philippe Devloo on 16/04/18.
//

#include "GeoNode.h"
#include "GeoElement.h"
#include <string>
#include "GeoMesh.h"
#include "tpanic.h"
#include "GeoElementSide.h"
#include "CompElementTemplate.h"
#include "GeoElementTemplate.h"
#include "CompMesh.h"
#include "Geom1d.h"
#include "GeomQuad.h"
#include "GeomTriangle.h"
#include "Shape1d.h"
#include "ShapeQuad.h"
#include "ShapeTriangle.h"
#include "DataTypes.h"


// Default constructor of CompMesh
CompMesh::CompMesh(){
    
}

// Default constructor of CompMesh
CompMesh::CompMesh(GeoMesh *geom){
    geomesh = geom;
}

// Copy constructor of CompMesh
CompMesh::CompMesh(const CompMesh &copy){
    
    geomesh = copy.geomesh;
    compelements=copy.compelements;
    dofs = copy.dofs;
    mathstatements = copy.mathstatements;
    solution = copy.solution;

}

// Destructor of CompMesh
CompMesh::~CompMesh(){
    
    
}

GeoMesh *CompMesh::GetGeoMesh() const{
    return geomesh;
}

void CompMesh::SetGeoMesh(GeoMesh *geo){
    geomesh = geo;
}

// Set the number of computational elements on the grid
void CompMesh::SetNumberElement(int64_t nelem){
    compelements.resize(nelem);
}

// Set the number of degrees of freedom
void CompMesh::SetNumberDOF(int64_t ndof){
    dofs.resize(ndof);
}

// Set the number of math statements
void CompMesh::SetNumberMath(int nmath){
    mathstatements.resize(nmath);
}

// Set the computational element associated to an index
void CompMesh::SetElement(int64_t elindex, CompElement *cel){
    compelements[elindex]=cel;
    
}
// Set the degree of freedom associated to an index
void CompMesh::SetDOF(int64_t index, const DOF &dof){
    dofs[index]=dof;
}

// Set the math statement object associated to an index
void CompMesh::SetMathStatement(int index, MathStatement *math){
    mathstatements[index]=math;
}

// Return the degree of freedom index
DOF &CompMesh::GetDOF(int64_t dofindex){
    return dofs[dofindex];
}

// Return the computational element associated to an index
CompElement *CompMesh::GetElement(int64_t elindex) const{
    return compelements[elindex];
}

// Return the math statement object associated to an index
MathStatement *CompMesh::GetMath(int matindex) const{
    return mathstatements[matindex];
}

// Return the vector with computational elements
std::vector<CompElement *> CompMesh::GetElementVec() const{
    return compelements;
}

// Return the vector with degrees of freedom
std::vector<DOF> CompMesh::GetDOFVec() const{
    return dofs;
}

// Return the vector with math statement objects
std::vector<MathStatement *> CompMesh::GetMathVec() const{
    return mathstatements;
}

// Set the vector with computational elements
void CompMesh::SetElementVec(const std::vector<CompElement *> &vec){
    compelements=vec;
}

// Set the vector with degrees of freedom
void CompMesh::SetDOFVec(const std::vector<DOF> &dofvec){
    dofs = dofvec;
}

// Set the vector with math statement objects
void CompMesh::SetMathVec(const std::vector<MathStatement *> &mathvec){
    mathstatements = mathvec;
}

// will create the computational elements
void CompMesh::AutoBuild(){
    int NumEl = geomesh->NumElements();
    for(int iel=0; iel <NumEl; iel++)
    {
    GeoElement *gel = GetGeoMesh()->Element(iel);
    CompElement *cel = gel->CreateCompEl(this, iel);
        SetNumberElement(iel+1);
        SetElement(iel, cel);
    }
    std::cout<<"Luego de autobuild debe de setarse el vectr de statements";
    
  
}

// Initialize the datastructure FirstEquation of the DOF objects
void CompMesh::Resequence(){
    std::cout<<"implementar Resequence en compmesh";
    DebugStop();
}

// Initialize the datastructure FirstEquation of the DOF objects in the order specified by the vector
void CompMesh::Resequence(VecInt &DOFindices){
    std::cout<<"implementar Resequence en compmesh";
    DebugStop();
}

std::vector<double> &CompMesh::Solution() {
    return solution;
}

void CompMesh::LoadSolution(std::vector<double> &Sol){
    solution = Sol;
}
