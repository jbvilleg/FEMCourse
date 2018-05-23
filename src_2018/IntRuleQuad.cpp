//
//  Geom1d.h
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.

#include "IntRuleQuad.h"
#include "IntRule1d.h"
//INT RULE ES UNA CLASE VIRTUAL
//Default Constructor of integration rule

IntRuleQuad::IntRuleQuad(){
    fOrder=0;
    SetOrder(fOrder);
}

IntRuleQuad::IntRuleQuad(int order){
    fOrder=order;
    SetOrder(order);
}

void IntRuleQuad::SetOrder(int order){
    fOrder =order;
    //numero de puntos necesarios para integrar con orden fOrder
    
    int npoints=0,resto=0;
    resto=fOrder%2;
    
    //ordem impar
    if (resto!=0) {
        npoints=(fOrder+1)/2;
    }
    if (resto==0) {
        npoints=(fOrder+2)/2;
    }
    
    npoints = npoints*npoints;
    
    fPoints.Resize(npoints,2);
    fWeights.resize(npoints);
    
    
    IntRule1d Int1Dx(fOrder);
    IntRule1d Int1Dy(fOrder);
    
    
    double weight;
    VecDouble co(2);
    for (int i=0; i<Int1Dx.NPoints(); i++) {
        
        Int1Dx.Point(i, co, weight);
        VecDouble coX(1);
        double weightX;
        coX[0]=co[0];
        weightX=weight;
        
        for (int j=0; j<Int1Dy.NPoints(); j++) {
            
            Int1Dy.Point(j, co, weight);
            
            fPoints(j+i*Int1Dy.NPoints(),0)=co[0];
            fPoints(j+i*Int1Dy.NPoints(),1)=coX[0];
            fWeights[j+i*Int1Dy.NPoints()]=weightX*weight;
        }
        
    }
    
    
}

void IntRuleQuad::gaulegQuad(const double x1, const double x2, VecDouble &x, VecDouble &w){
    
    IntRule1d IntGauss1Dx(fOrder);
    IntRule1d IntGauss1Dy(fOrder);
    double nPoints = x.size();
    VecDouble weightx(x.size()), coX(nPoints);
    VecDouble weighty(x.size()), coY(nPoints);
    
    IntGauss1Dx.gauleg(x1, x2, coX, weightx);
    IntGauss1Dy.gauleg(x1, x2, coY, weighty);
    
    x.resize(2*nPoints*nPoints);
    w.resize(nPoints*nPoints);
    
    for (int i = 0; i<nPoints; i++) {
        
        for (int j = 0; j<nPoints; j++) {
            w[j+i*nPoints]=weightx[j]*weighty[i];
            x[j+i*nPoints]=coX[j];
            x[j+i*nPoints+nPoints*nPoints]=coY[i];
        }
    }
    
    
    
}
