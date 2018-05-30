//
//  Poisson.h
//  FemCourse
//
//  Created by Philippe Devloo on 24/04/18.
//

#include "Poisson.h"
#include "MathStatement.h"
#include "tpanic.h"

Poisson::Poisson(){
    
}

Poisson::Poisson(Matrix &perm){
    permeability=perm;
}

Poisson::Poisson(const Poisson &copy){
    permeability=copy.permeability;
    forceFunction=copy.forceFunction;
}

Poisson &Poisson::operator=(const Poisson &copy){
    permeability=copy.permeability;
    forceFunction=copy.forceFunction;
    return *this;
}

Poisson *Poisson::Clone() const{
   
}

Poisson::~Poisson(){
    
}

Matrix Poisson::GetPermeability() const{
    return permeability;
}

void Poisson::SetPermeability(const Matrix &perm){
   
    permeability=perm;
}

int Poisson::NState() const{
    return 1;
}

// Method to implement integral over element's volume
void Poisson::Contribute(IntPointData &data, double weight , Matrix &EK, Matrix &EF) const{
    weight = data.weight;
    
    
    
}

void Poisson::ContributeError(IntPointData &integrationpointdata, std::function<void (const VecDouble &, VecDouble &, Matrix &)> &exact){
    
}
// Prepare and print post processing data
VecDouble Poisson::PostProcess(const IntPointData &integrationpointdata, const PostProcVar var) const{


}

