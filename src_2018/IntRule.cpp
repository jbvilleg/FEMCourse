//
//  Geom1d.h
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//


#include "IntRule.h"
#include "IntRule1d.h"
//INT RULE ES UNA CLASE VIRTUAL
//Default Constructor of integration rule
IntRule::IntRule(){
    fOrder = 0;
    SetOrder(fOrder);
}

//Constructor of integration rule
IntRule::IntRule(int order){
    fOrder=order;
    SetOrder(fOrder);
}
//Destructor of integration rule
IntRule::~IntRule(){
    
}

//Copy constructor of integration rule
IntRule::IntRule(const IntRule &copy){

    fWeights=copy.fWeights;
    fPoints = copy.fPoints;
    fOrder=copy.fOrder;
    SetOrder(fOrder);
    
}



//operator=
IntRule &IntRule::operator=(const IntRule &copy){
    fPoints=copy.fPoints;
    fWeights=copy.fWeights;
    fOrder=copy.fOrder;
    SetOrder(fOrder);
}


//Method to return the number of integration points
int IntRule::NPoints() const{
    //el numero de puntos es calculado por el tamaño del vector
    int npoints = fWeights.size();
    return npoints;
}

//Function returning coordinates and weights of integration points
void IntRule::Point(int p, VecDouble &co, double &weight) const {
    co.resize(2);
    co[0]=fPoints.GetVal(p, 0);
    co[1]=fPoints.GetVal(p, 1);
    weight=fWeights[p];
    
}





//Fuction for printing results
void IntRule::Print(std::ostream &out) const{
    std::cout<<"Orden Polinomial : "<<fOrder<<"\n";
    std::cout<<"Numero de Puntos de Integración: "<<NPoints()<<"\n";

}
