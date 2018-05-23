
#include "MathStatement.h"
#include "DataTypes.h"
#include  "IntPointData.h"
#include "Poisson.h"
#include "DOF.h"


// Default constructor of DOF
DOF::DOF(){
    
}

// Copy constructor of DOF
DOF::DOF(const DOF &copy){
    // First associated equation
    
    firstequation = copy.firstequation;
    nshape = copy.nshape;
    nstate = copy.nstate;
    order = copy.order;
}

// Operator of copy
DOF &DOF::operator=(const DOF &copy){
    firstequation = copy.firstequation;
    nshape = copy.nshape;
    nstate = copy.nstate;
    order = copy.order;
    return *this;
}

// Destructor of DOF
DOF::~DOF(){
    
}

// Return the first associated equation
int64_t DOF::GetFirstEquation(){
    return firstequation;
}

// Set the first associated equation
void DOF::SetFirstEquation(int64_t first){
    
    firstequation=first;

}

// Set the number of shape functions associated with the DOF and state variables associated with each shape function
void DOF::SetNShapeStateOrder(int NShape, int NState, int Order){
    nshape = NShape;
    nstate = NState;
    order = Order;
}

// Return the number of shape functions associated with the DOF
int DOF::GetNShape() const{
    return nshape;
}

// Return the number of state variables associated with each shape function
int DOF::GetNState() const{
    return nstate;
}

// Return maximum the order of the polinoms associated with DOF
int DOF::GetOrder() const{
    return order;
}

