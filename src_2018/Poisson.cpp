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
    // return new MathStatement(*this);
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
   
    VecDouble &phi=data.phi;
    Matrix &dphi=data.dphidx;
    VecDouble &x = data.x;
    Matrix &axes =data.axes;
    
    int nphi= phi.size();
    int nphij=phi.size();
    VecDouble f(nphi,0.);
    forceFunction(x,f);
    Matrix perm = GetPermeability();
    Matrix Kdphi(2,2,0.);
    
    
    for (int i = 0; i<nphi; i++) {
        EF(i,0)+= -weight*f[i]*phi[i];
        for(int j = 0; j<nphi; j++){
            //revisar esto! k¿?
            for (int k=0; k<nphi; k++) {
                EK(i,j) += weight*(dphi(k,i)*dphi(k,j));
            }
            
        }
    }
    
    
}

// Prepare and print post processing data
VecDouble Poisson::PostProcess(const IntPointData &integrationpointdata, const PostProcVar var) const{


}

