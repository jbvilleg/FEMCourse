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
    
    for(int i =0; i<=phi.size(); i++){
        x[0] += phi[i]*NodeCo(0,0); //verificar la matriz;
        x[1] += phi[i]*NodeCo(0,1); //verificar la matriz;
        x[2] += 0;
    }
 
}

/// Computes the value of x and gradx for a given point in parameter space
void GeomQuad::GradX(const VecDouble &xi, Matrix &NodeCo, VecDouble &x, Matrix &gradx){
    
    VecDouble phi(4);
    Matrix dphi(2,3);
    x.resize(3);
    Shape(xi, phi, dphi);
    X(xi,NodeCo,x);
    
    for(int i=0; i<= 3; i++){
        for(int j=0; j<NodeCo.Rows();j++){
            gradx(i,0)+= NodeCo(j,i)*dphi(0,i);
            gradx(i,1)+= NodeCo(j,i)*dphi(1,i);
        }
    }
    
    
}

/// return the number of nodes of the template
int GeomQuad::NumNodes(){
    //revisar!
    //como determinar el numero de nodos¿?
    //usar topology?
    
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
void GeomQuad::SetNeighbour(int side, GeoElementSide &neighbour){
    fNeighbours[side] = neighbour;
}
