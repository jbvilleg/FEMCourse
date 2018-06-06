//
//  Poisson.h
//  FemCourse
//
//  Created by Philippe Devloo on 24/04/18.
//

#include "Poisson.h"
#include "MathStatement.h"
#include "tpanic.h"
#include <functional>

Poisson::Poisson(){
    
}

Poisson::Poisson(int materialid, Matrix &perm){
    SetMatID(materialid);
    permeability=perm;
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
    return 2;
}



// Method to implement integral over element's volume
void Poisson::Contribute(IntPointData &data,double weight ,Matrix &EK, Matrix &EF) const{
    
    VecDouble phi = data.phi;
    Matrix dphi = data.dphidx;
    Matrix axes = data.axes;
    
    int nshape = phi.size();
    int dim = dphi.Rows();
    
    std::vector<double> du(dim);
    std::vector<double> dv(dim);

   Matrix perm(dim, dim);
   perm = permeability;

    for (int i = 0; i < nshape; i++) {
        dv[0] = dphi(0, i) * axes(0, 0) + dphi(1, i) * axes(1, 0);
        dv[1] = dphi(0, i) * axes(0, 1) + dphi(1, i) * axes(1, 1);
        
        VecDouble f(2,0.);
        forceFunction(data.x,f);
         //se n statate var for 1 (u)
//        EF(i,0)+= weight*phi[i]*f[0];
//        EF(i,0)+= weight*phi[i]*f[1];

        //se n statate var for 2 (ux uy)
        EF(2*i,0)+= weight*phi[i]*f[0];
        EF(2*i+1,0)+= weight*phi[i]*f[1];
        
        for (int j = 0; j < nshape; j++) {
            du[0] = dphi(0, j) * axes(0, 0) + dphi(1, j) * axes(1, 0);
            du[1] = dphi(0, j) * axes(0, 1) + dphi(1, j) * axes(1, 1);

            //se n statate var for 1 (u)
//            EK(i, j) += weight*(du[0] * dv[0] * perm(0, 0)  + du[1] * dv[0] * perm(0, 1));
//            EK(i, j) += weight*(du[0] * dv[1] * perm(1, 0)  + du[1] * dv[1] * perm(1, 1));
            
            //se n statate var for 2 (ux uy)
            EK(2*i,2*j) += weight*(du[0]*dv[0]*perm(0,0)+du[0]*dv[1]*perm(1,0));
            EK(2*i,2*j+1) += weight*(du[0]*dv[0]*perm(0,1)+du[0]*dv[1]*perm(1,1));
            EK(2*i+1,2*j) += weight*(du[1]*dv[0]*perm(0,0)+du[1]*dv[1]*perm(1,0));
            EK(2*i+1,2*j+1) += weight*(du[1]*dv[0]*perm(0,1)+du[1]*dv[1]*perm(1,1));

        }
     
    }
}

void Poisson::ContributeError(IntPointData &integrationpointdata, std::function<void (const VecDouble &, VecDouble &, Matrix &)> &exact){
    
}
// Prepare and print post processing data
VecDouble Poisson::PostProcess(const IntPointData &integrationpointdata, const PostProcVar var) const{


}

