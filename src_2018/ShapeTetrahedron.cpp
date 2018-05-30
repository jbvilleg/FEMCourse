    //
//  Geom1d.h
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.


#include "IntRuleTetrahedron.h"
#include "ShapeTetrahedron.h"
#include "TopologyTetrahedron.h"
#include "tpanic.h"

/// computes the shape functions in function of the coordinate in parameter space and orders of the shape functions (size of orders is number of sides of the element topology)

 void ShapeTetrahedron::Shape(const VecDouble &xi, VecInt &orders, VecDouble &phi, Matrix &dphi){
    
     int nshape = NShapeFunctions(orders);
     //int nsides = orders.size();
     
     phi.resize(nshape);
     dphi.Resize(2, nshape);
     
     if (nshape==4) {
         phi[0]=1-xi[0]-xi[1]-xi[2];
         phi[1]=xi[0];
         phi[2]=xi[1];
         phi[3]=xi[2];
         
         dphi(0,0)=-1;
         dphi(0,1)=1;
         dphi(0,2)=0;
         dphi(0,3)=0;
         
         dphi(1,0)=-1;
         dphi(1,1)=0;
         dphi(1,2)=1;
         dphi(1,2)=0;
         
         dphi(1,0)=-1;
         dphi(1,1)=0;
         dphi(1,2)=0;
         dphi(1,2)=1;
     }
     
     if (nshape==10) {
         
         VecDouble eps(3);
         eps[0]=1-xi[0]-xi[1]-xi[2];
         eps[1]=xi[0];
         eps[2]=xi[1];
         eps[3]=xi[2];
         
         
         for (int i=0; i<=3; i++) {
             phi[i]=eps[i]*(2.*eps[i]-1.);
         }
         phi[4]=4.*phi[0]*phi[1];
         phi[5]=4.*phi[1]*phi[2];
         phi[6]=4.*phi[2]*phi[0];
         phi[7]=4.*phi[0]*phi[3];
         phi[8]=4.*phi[1]*phi[3];
         phi[9]=4.*phi[2]*phi[3];
         
         dphi(0,0)=1.-4.*(1-eps[0]-eps[1]-eps[2]);
         dphi(0,1)=-1.+4.*eps[0];
         dphi(0,2)=0.;
         dphi(0,3)=0.;
         dphi(0,4)= -4.*eps[0]*(-1.+2.*eps[0])*(-1.+2.*(1.-eps[0]-eps[1]-eps[2]))-8.*eps[0]*(-1.+2.*eps[0])*(1.-eps[0]-eps[1]-eps[2])+8.*eps[0]*(-1.+2.*(1.-eps[0]-eps[1]-eps[2]))*(1.-eps[0]-eps[1]-eps[2])+4.*(-1.+2.*eps[0])*(-1.+2.*(1.-eps[0]-eps[1]-eps[2]))*(1.-eps[0]-eps[1]-eps[2]);
         dphi(0,5)= 8.*eps[0]*eps[1]*(-1.+2.*eps[1])+4.*(-1.+2.*eps[0])*eps[1]*(-1.+2.*eps[1]);
         dphi(0,6)= -4*eps[1]*(-1.+2.*eps[1])*(-1.+2.*(1.-eps[0]-eps[1]-eps[2]))-8.*eps[1]*(-1.+2.*eps[1])*(1.-eps[0]-eps[1]-eps[2]);
         dphi(0,7)= -4*(-1.+2.*(1.-eps[0]-eps[1]-eps[2]))*eps[2]*(-1.+2.*eps[2])-8.*(1.-eps[0]-eps[1]-eps[2])*eps[2]*(-1.+2.*eps[2]);
         dphi(0,8)= 8.*eps[0]*eps[2]*(-1.+2.*eps[2])+4.*(-1.+2.*eps[0])*eps[2]*(-1.+2.*eps[2]);
         dphi(0,9)=0.;
         
         dphi(1,0)= 1.-4.*(1.-eps[0]-eps[1]-eps[2]);
         dphi(1,1)= 0.;
         dphi(1,2)= -1.+4.*eps[1];
         dphi(1,3)= 0.;
         dphi(1,4)= -4.*eps[0]*(-1.+2.*eps[0])*(-1.+2.*(1.-eps[0]-eps[1]-eps[2]))-8.*eps[0]*(-1.+2.*eps[0])*(1.-eps[0]-eps[1]-eps[2]);
         dphi(1,5)= 8.*eps[0]*(-1.+2.*eps[0])*eps[1]+4*eps[0]*(-1.+2.*eps[0])*(-1+2.*eps[1]);
         dphi(1,6)= -4.*eps[1]*(-1.+2.*eps[1])*(-1.+2.*(1.-eps[0]-eps[1]-eps[2]))-8.*eps[1]*(-1.+2.*eps[1])*(1.-eps[0]-eps[1]-eps[2])+8.*eps[1]*(-1.+2.*(1.-eps[0]-eps[1]-eps[2]))*(1.-eps[0]-eps[1]-eps[2])+4.*(-1.+2.*eps[1])*(-1.+2.*(1.-eps[0]-eps[1]-eps[2]))*(1.-eps[0]-eps[1]-eps[2]);
         dphi(1,7)=-4.*(-1.+2.*(1.-eps[0]-eps[1]-eps[2]))*eps[2]*(-1.+2.*eps[2])-8.*(1.-eps[0]-eps[1]-eps[2])*eps[2]*(-1.+2.*eps[2]);
         dphi(1,8)=0.;
         dphi(1,9)= 8.*eps[1]*eps[2]*(-1.+2.*eps[2])+4.*(-1.+2.*eps[1])*eps[2]*(-1.+2.*eps[2]);
         
         dphi(2,0)= 1.-4.*(1.-eps[0]-eps[1]-eps[2]);
         dphi(2,1)= 0.;
         dphi(2,2)= 0.;
         dphi(2,3)= -1.+4.*eps[2];
         dphi(2,4)= -4.*eps[0]*(-1.+2.*eps[0])*(-1.+2.*(1.-eps[0]-eps[1]-eps[2]))-8.*eps[0]*(-1.+2.*eps[0])*(1.-eps[0]-eps[1]-eps[2]);
         dphi(2,5)= 0.;
         dphi(2,6)= -4*eps[1]*(-1.+2.*eps[1])*(-1.+2.*(1.-eps[0]-eps[1]-eps[2]))-8.*eps[1]*(-1.+2.*eps[1])*(1.-eps[0]-eps[1]-eps[2]);
         dphi(2,7)= 8.*(-1.+2.*(1.-eps[0]-eps[1]-eps[2]))*(1.-eps[0]-eps[1]-eps[2])*eps[2]+4.*(-1.+2.*(1.-eps[0]-eps[1]-eps[2]))*(1.-eps[0]-eps[1]-eps[2])*(-1.+2.*eps[2])-4.*(-1+2.*(1.-eps[0]-eps[1]-eps[2]))*eps[2]*(-1.+2.*eps[2])-8.*(1.-eps[0]-eps[1]-eps[2])*eps[2]*(-1.+2.*eps[2]);
         dphi(2,8)= 8.*eps[0]*(-1.+2.*eps[0])*eps[2]+4.*eps[0]*(-1.+2.*eps[0])*(-1.+2.*eps[2]);
         dphi(2,9)= 8.*eps[1]*(-1.+2.*eps[1])*eps[2]+4.*eps[1]*(-1.+2.*eps[1])*(-1.+2.*eps[2]);
     }
     
     
     
}

/// returns the number of shape functions associated with a side
 int ShapeTetrahedron::NShapeFunctions(int side, int orders){
     if (side<4) {
         return 1;
     }
     else{
         return orders-1;
     }
}

/// returns the total number of shape functions
 int ShapeTetrahedron::NShapeFunctions(VecInt &orders){
     int nSides = orders.size();
     int nshape=0;
     for(int i=0; i<= nSides-1;i++){
         nshape = nshape + NShapeFunctions(i,orders[i]);
     }
     return nshape;
}
