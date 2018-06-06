//
//  Geom1d.h
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//


#include "Topology1d.h"
#include "Geom1d.h"
#include "GeomQuad.h"
#include "GeoElementSide.h"

/// Constructor
GeomQuad::GeomQuad(){
    fNodeIndices.resize(nCorners);
    for(int i=0; i<nCorners; i++){
        fNodeIndices[i]=-1;
    }
}

/// destructor
GeomQuad::~GeomQuad(){
    
}

/// copy constructor
GeomQuad::GeomQuad(const GeomQuad &copy){
    fNodeIndices=copy.fNodeIndices;
    for(int i =0; i<=nSides; i++){
        fNeighbours[i]=copy.fNeighbours[i];
    }
    
}

/// operator=
GeomQuad &GeomQuad::operator=(const GeomQuad &copy){
    
    fNodeIndices=copy.fNodeIndices;
    for(int i =0; i<=nSides; i++){
        fNeighbours[i]=copy.fNeighbours[i];
    }
    return *this;
}

/// Computes the shape functions associated with the geometric map
void GeomQuad::Shape(const VecDouble &xi, VecDouble &phi, Matrix &dphi){
   
    phi[0] = 0.25*(1.-xi[0])*(1.-xi[1]);
    phi[1] = 0.25*(1.+xi[0])*(1.-xi[1]);
    phi[2] = 0.25*(1.+xi[0])*(1.+xi[1]);
    phi[3] = 0.25*(1.-xi[0])*(1.+xi[1]);
    
    dphi(0,0) = 0.25*(xi[1]-1.);
    dphi(1,0) = 0.25*(xi[0]-1.);
    
    dphi(0,1) = 0.25*(1.-xi[1]);
    dphi(1,1) =-0.25*(1.+xi[0]);
    
    dphi(0,2) = 0.25*(1.+xi[1]);
    dphi(1,2) = 0.25*(1.+xi[0]);
    
    dphi(0,3) =-0.25*(1.+xi[1]);
    dphi(1,3) = 0.25*(1.-xi[0]);
}

/// Computes the value of x for a given point in parameter space as a function of corner coordinates
void GeomQuad::X(const VecDouble &xi, Matrix &NodeCo, VecDouble &x){
    
    VecDouble phi(4);
    Matrix dphi(2,4);
    x.resize(3);
    Shape(xi, phi, dphi);
    
    int nSpace = NodeCo.Rows();
    for(int i = 0; i < nSpace; i++) {
        x[i] = 0.0;
        for(int j = 0; j < 4; j++) {
            x[i] += phi[j]*NodeCo.GetVal(i,j);
        }
    }
    
//    for(int i =0; i<=phi.size(); i++){
//        x[0] += phi[i]*NodeCo.GetVal(0,0); //verificar la matriz;
//        x[1] += phi[i]*NodeCo.GetVal(0,1); //verificar la matriz;
//        x[2] += 0;
//    }

}

/// Computes the value of x and gradx for a given point in parameter space
void GeomQuad::GradX(const VecDouble &xi, Matrix &NodeCo, VecDouble &x, Matrix &gradx){
    
    VecDouble phi(4,0.0);
    Matrix dphi(2,4,0.0);
    int dim = Dimension;
    Shape(xi, phi, dphi);
    gradx.Resize(dim,2);
    gradx.Zero();
    
    for(int i=0; i< 4; i++){
        for(int j=0; j<Dimension;j++){
            gradx(j,0)+= NodeCo(j,i)*dphi(0,i);
            gradx(j,1)+= NodeCo(j,i)*dphi(1,i);
        }
    }
    
    
}

/// return the number of nodes of the template
int GeomQuad::NumNodes(){
    
    return nCorners;
}

/// Set the node indices of the element
void GeomQuad::SetNodes(const VecInt &nodes){
    fNodeIndices =nodes;
}

/// Set the node indices of the element
void GeomQuad::GetNodes(VecInt &nodes){
    nodes = fNodeIndices;
}

/// Return the index of a node
int GeomQuad::NodeIndex(int node){
    return fNodeIndices[node];
}

/// Return the neighbour along side
GeoElementSide GeomQuad::Neighbour(int side){
    return fNeighbours[side];
}

/// Initialize the neighbour data structure
void GeomQuad::SetNeighbour(int side, const GeoElementSide &neighbour){
    fNeighbours[side] = neighbour;
}
