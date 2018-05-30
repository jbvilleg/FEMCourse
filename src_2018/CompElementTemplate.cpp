
#include "CompElement.h"
#include "CompElementTemplate.h"
#include "CompMesh.h"

// Default constructor of CompElementTemplate
template<class Shape>
CompElementTemplate<Shape>::CompElementTemplate():CompElement(){
    std::cout<<"bingo";
    
    
    
}

template<class Shape>
CompElementTemplate<Shape>::CompElementTemplate(int64_t ind, CompMesh *cmesh,  GeoElement *geo) : CompElement(ind,cmesh,geo){
    
    //numero de elementos en la malla computacional
    int Nelem = cmesh->GetElementVec().size();
    cmesh->SetNumberElement(Nelem+1);
    //recordar que inicia la numeracion en 0
    cmesh->SetElement(Nelem, this);
    this->SetIndex(ind);
    intrule.SetOrder(1);
   
}

// Copy constructor of CompElementTemplate
template<class Shape>
CompElementTemplate<Shape>::CompElementTemplate(const CompElementTemplate &){
    
}

// Operator of copy
template<class Shape>
CompElementTemplate<Shape> &CompElementTemplate<Shape>::operator=(const CompElementTemplate &){
    
    
}

// Destructor of CompElementTemplate
template<class Shape>
CompElementTemplate<Shape>::~CompElementTemplate(){
    
}

// Method for creating a copy of the element
template<class Shape>
CompElement * CompElementTemplate<Shape>::Clone() const{
    
    CompElement *cel;
    return cel;
}

// Compute shape functions set at point x
template<class Shape>
void CompElementTemplate<Shape>::ShapeFunctions(const VecDouble &intpoint, VecDouble &phi, Matrix &dphi) const {
    
}

// Return the number of shape functions
template<class Shape>
int CompElementTemplate<Shape>::NShapeFunctions() const{
    int order = intrule.GetOrder();
    VecInt orders(1, order);
    return Shape::NShapeFunctions(orders);
}
// Return the number of degree of freedom
template<class Shape>
int CompElementTemplate<Shape>::NDOF() const {
    
    return -1;
}

// Return the number of shape functions stored in the DOF data structure
template<class Shape>
int CompElementTemplate<Shape>::NShapeFunctions(int doflocindex){
    
    return -1;
}

// Use the Shape template class to compute the number of shape functions
template<class Shape>
int CompElementTemplate<Shape>::ComputeNShapeFunctions(int doflocindex, int order){
    
    return -1;
}

// Return space dimension
template<class Shape>
int CompElementTemplate<Shape>::Dimension() const {
    std::cout<<"revisar Dimension en Shape";
    return 1;
}

template class CompElementTemplate<Shape1d>;
template class CompElementTemplate<ShapeQuad>;
template class CompElementTemplate<ShapeTriangle>;
template class CompElementTemplate<ShapeTetrahedron>;
