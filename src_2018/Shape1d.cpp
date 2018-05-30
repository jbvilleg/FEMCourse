    //
//  Geom1d.h
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.

#include "Shape1d.h"
#include "tpanic.h"

/// computes the shape functions in function of the coordinate in parameter space and orders of the shape functions (size of orders is number of sides of the element topology)

void Shape1d::Shape(const VecDouble &xi, VecInt &orders, VecDouble &phi, Matrix &dphi){

    
    int nshape = NShapeFunctions(orders);
    int order = nshape - 1;
    
    phi.resize(nshape);
    dphi.Resize(1, nshape);
    
    for (int i=0; i<nshape; i++) {
        phi[i]=1.;
        dphi(0,i)=0;
    }
    
    for (int i=0; i<nshape; i++) {
        double epsi=-1.+i*2./order;
        
        for (int j=0; j<nshape; j++) {
            
            Matrix axdphi(1,nshape);
            if (i!=j) {
                double epsj=-1.+j*2./order;
                
                phi[i]*=(xi[0]-epsj)/(epsi-epsj);
                
                axdphi(0,i)=1/(epsi-epsj);
                
                for (int k=0; k<nshape; k++) {
                    if (k!=i&&k!=j) {
                        epsj=-1.+k*2./order;
                        axdphi(0,i)*=(xi[0]-epsj)/(epsi-epsj);
                    }
                    
                }
                
                dphi(0,i)+=axdphi(0,i);
                
            }
            
        }
        
    }
    
    
}

/// returns the number of shape functions associated with a side

int Shape1d::NShapeFunctions(int side, int orders){
    
    if (side<2) {
        return 1;
    }
    else{
        return orders-1;
    }

}

/// returns the total number of shape functions
int Shape1d::NShapeFunctions(VecInt &orders){
    int nSides = orders.size();
    int nshape=0;
    for(int i=0;i<=nSides-1;i++){
        nshape = nshape + NShapeFunctions(i,orders[i]);
    }
    return nshape;
    
}
