//
//  Geom1d.h
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//


#include "Topology1d.h"
#include "Geom1d.h"
#include "GeomQuad.h"
#include "GeomTriangle.h"
#include "GeoElementSide.h"

/// Constructor
GeomTriangle::GeomTriangle(){
    
}

/// destructor
GeomTriangle::~GeomTriangle(){
    
}

/// copy constructor
GeomTriangle::GeomTriangle(const GeomTriangle &copy){
    fNodeIndices = copy.fNodeIndices;
    for(int i =0; i<=nSides; i++){
        fNeighbours[i]=copy.fNeighbours[i];
    }
}

/// operator=
GeomTriangle &GeomTriangle::operator=(const GeomTriangle &copy){
    fNodeIndices = copy.fNodeIndices;
    for(int i =0; i<=nSides; i++){
        fNeighbours[i]=copy.fNeighbours[i];
    }
}

/// Computes the shape functions associated with the geometric map
 void GeomTriangle::Shape(const VecDouble &xi, VecDouble &phi, Matrix &dphi){
    
      //xi = xi[0]
      //eta = xi[1];
     phi.resize(3);
     dphi.Resize(2, 3);
     phi[0] = 1.0-xi[0]-xi[1];
     phi[1] = xi[0];
     phi[2] = xi[1];
     
     dphi(0,0) = dphi(1,0) = -1.0;
     dphi(0,1) = dphi(1,2) =  1.0;
     dphi(1,1) = dphi(0,2) =  0.0;
     
}

/// Computes the value of x for a given point in parameter space as a function of corner coordinates
 void GeomTriangle::X(const VecDouble &xi, Matrix &NodeCo, VecDouble &x){
     VecDouble phi(3);
     Matrix dphi(2,3);
     x.resize(3);
     
     Shape(xi, phi, dphi);
     for(int i =0; i<=phi.size(); i++){
         x[0] += phi[i]*NodeCo(0,0); //verificar la matriz;
         x[1] += phi[i]*NodeCo(0,1); //verificar la matriz;
         x[2] += 0;
     }
}

/// Computes the value of x and gradx for a given point in parameter space
 void GeomTriangle::GradX(const VecDouble &xi, Matrix &NodeCo, VecDouble &x, Matrix &gradx){
     
     VecDouble phi(3,0);
     Matrix dphi(2,3);
     X(xi,NodeCo,x);
     Shape(xi,phi,dphi);
     
     for(int i=0; i< 3; i++){
         for(int j=0; j<NodeCo.Rows();j++){
         gradx(i,0)+= NodeCo(j,i)*dphi(0,i);
         gradx(i,1)+= NodeCo(j,i)*dphi(1,i);
         }
     }
     
     
}

/// return the number of nodes of the template
int GeomTriangle::NumNodes(){
    //verificar
    return nCorners;
}

/// Set the node indices of the element
void GeomTriangle::SetNodes(const VecInt &nodes){
    fNodeIndices = nodes;
}

/// Set the node indices of the element
void GeomTriangle::GetNodes(VecInt &nodes){
    nodes = fNodeIndices;
}

/// Return the index of a node
int GeomTriangle::NodeIndex(int node){
    node = fNodeIndices[node];
}

/// Return the neighbour along side
GeoElementSide GeomTriangle::Neighbour(int side){
    return fNeighbours[side];
}

/// Initialize the neighbour data structure
void GeomTriangle::SetNeighbour(int side, GeoElementSide &neighbour){
    fNeighbours[side] = neighbour;
}
