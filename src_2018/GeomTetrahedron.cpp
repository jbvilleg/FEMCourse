//
//  GeomTetrahedron.cpp
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//

#include "TopologyTetrahedron.h"
#include "GeomTetrahedron.h"

    /// Constructor
    GeomTetrahedron::GeomTetrahedron(){
        fNodeIndices.resize(nCorners);
        for(int i=0; i<nCorners; i++) fNodeIndices[i]=-1;
    }
    
    /// destructor
    GeomTetrahedron::~GeomTetrahedron(){
        
    }
    
    /// copy constructor
    GeomTetrahedron::GeomTetrahedron(const GeomTetrahedron &copy){
        fNodeIndices=copy.fNodeIndices;
        
        for (int i=0; i<nSides; i++) {
            fNeighbours[i]=copy.fNeighbours[i];
        }
    }
    
    /// operator=
    GeomTetrahedron &GeomTetrahedron::operator=(const GeomTetrahedron &copy){
        fNodeIndices=copy.fNodeIndices;
        
        for (int i=0; i<nSides; i++) {
            fNeighbours[i]=copy.fNeighbours[i];
        }
        return *this;
    }
    
    /// Computes the shape functions associated with the geometric map
    void GeomTetrahedron::Shape(const VecDouble &xi, VecDouble &phi, Matrix &dphi){
        double qsi = xi[0], eta = xi[1], zeta  = xi[2];
        phi[0]  = 1.0-qsi-eta-zeta;
        phi[1]  = qsi;
        phi[2]  = eta;
        phi[3]  = zeta;
        
        dphi(0,0) = -1.0;
        dphi(1,0) = -1.0;
        dphi(2,0) = -1.0;
        dphi(0,1) =  1.0;
        dphi(1,1) =  0.0;
        dphi(2,1) =  0.0;
        dphi(0,2) =  0.0;
        dphi(1,2) =  1.0;
        dphi(2,2) =  0.0;
        dphi(0,3) =  0.0;
        dphi(1,3) =  0.0;
        dphi(2,3) =  1.0;
    }
    
    /// Computes the value of x for a given point in parameter space as a function of corner coordinates
    void GeomTetrahedron::X(const VecDouble &xi, Matrix &NodeCo, VecDouble &x){
        
        VecDouble phi(4,0.);
        Matrix dphi(3,4,0.);
        Shape(xi,phi,dphi);
        int space = NodeCo.Rows();
        for(int i = 0; i < space; i++) {
            x[i] = 0.0;
            for(int j = 0; j < 4; j++) {
                x[i] += phi[j]*NodeCo(i,j);
            }
        }
        
    }
    
    /// Computes the value of x and gradx for a given point in parameter space
    void GeomTetrahedron::GradX(const VecDouble &xi, Matrix &NodeCo, VecDouble &x, Matrix &gradx){
        
        gradx.Resize(3,3);
        gradx.Zero();


    // revisar Geo Tetrahedron
        
        VecDouble phi(4,0.);
        Matrix dphi(3,4);
        X(xi,NodeCo,x);
        Shape(xi,phi,dphi);
        for(int i = 0; i < 4; i++)
        {
            for(int j = 0; j < 3; j++)
            {
                gradx(j,0) += NodeCo(j,i)*dphi(0,i);
                gradx(j,1) += NodeCo(j,i)*dphi(1,i);
                gradx(j,2) += NodeCo(j,i)*dphi(2,i);
            }
        }
        
    }

    /// return the number of nodes of the template
    int GeomTetrahedron::NumNodes(){
        return nCorners;
    }

    /// Set the node indices of the element
    void GeomTetrahedron::SetNodes(const VecInt &nodes){
        fNodeIndices=nodes;
    }
    
    /// Set the node indices of the element
    void GeomTetrahedron::GetNodes(VecInt &nodes){
        nodes=fNodeIndices;
    }
    
    /// Return the index of a node
    int GeomTetrahedron::NodeIndex(int node){
        return fNodeIndices[node];
    }

    /// Return the neighbour along side
    GeoElementSide GeomTetrahedron::Neighbour(int side){
        return fNeighbours[side];
    }

    /// Initialize the neighbour data structure
    void GeomTetrahedron::SetNeighbour(int side, const GeoElementSide &neighbour){
        fNeighbours[side]=neighbour;
    }
