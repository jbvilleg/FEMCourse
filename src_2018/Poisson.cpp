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

int Poisson::NEvalErrors() const{
    return 3;
}

int Poisson::NState() const{
    return 1;
}

int Poisson::VariableIndex(const PostProcVar var) const{
    
    int nvar = 10;
    for (int i=0; i<nvar; i++) {
        if (var==PostProcVar(i)) {
            return i;
        }
    }
    return 0;
}

// Return the variable index associated with the name
Poisson::PostProcVar Poisson::VariableIndex(const std::string &name){
    
    if (!strcmp("Sol", name.c_str()))  return ESol;
    if (!strcmp("Solution", name.c_str()))  return ESol;
    if (!strcmp("DSol", name.c_str()))  return EDSol;
    if (!strcmp("DSolution", name.c_str()))  return EDSol;
    if (!strcmp("Flux", name.c_str()))         return EFlux;
    if (!strcmp("Force", name.c_str()))   return EForce;
    if (!strcmp("Sol_exact", name.c_str()))   return ESolExact;
    if (!strcmp("DSol_exact", name.c_str()))   return EDSolExact;
    
    std::cout  << " Var index not implemented " << std::endl;
    DebugStop();
    return ENone;
}

// Return the number of variables associated with the variable indexed by var. Param var Index variable into the solution, is obtained by calling VariableIndex
int Poisson::NSolutionVariables(const PostProcVar var){
    
    switch(var) {
            
        case ESol:
            return this->Dimension(); // Solution, Vector
        case EDSol:
            return this->Dimension();// Derivative of solution, Vector
        case EFlux:
            return this->Dimension(); // Flux, Vector
        case EForce:
            return this->Dimension(); // Force vector, Vector
        case ESolExact:
            return this->Dimension(); // Sol_exact, Vector
        case EDSolExact:
            return this->Dimension(); // DSol_exact, Vector
            
        default:
        {
            std::cout  << " Var index not implemented " << std::endl;
            DebugStop();
        }
    }
    return 0;
    
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
     //    se n statate var for 1 (u)
        EF(i,0)+= weight*phi[i]*f[0];

        
        for (int j = 0; j < nshape; j++) {
            du[0] = dphi(0, j) * axes(0, 0) + dphi(1, j) * axes(1, 0);
            du[1] = dphi(0, j) * axes(0, 1) + dphi(1, j) * axes(1, 1);

            //se n statate var for 1 (u)
            EK(i, j) += weight*(du[0] * dv[0] * perm(0, 0)  + du[1] * dv[0] * perm(0, 1));
            EK(i, j) += weight*(du[0] * dv[1] * perm(1, 0)  + du[1] * dv[1] * perm(1, 1));
        }
     
    }
}

void Poisson::ContributeError(IntPointData &integrationpointdata, VecDouble &u_exact, Matrix &du_exact, VecDouble &errors) const{
   
    
    
    
    
    
    
}

std::vector<double> Poisson::PostProcessSolution(const IntPointData &data, const int varindex) const{
    
    PostProcVar var = PostProcVar(varindex);
    
    VecDouble u_h = data.solution;
  
    Matrix du_h = data.dsoldx;
   
    Matrix perm = GetPermeability();
    
    VecDouble Solout;
    switch(var) {
            
        case ESol: //u
        {
            
            Solout.resize(2);
            Solout[0] = u_h[0];
          
            
        }
            break;
            
       
            
        default:
        {
            std::cout  << " Var index not implemented " << std::endl;
            DebugStop();
        }
    }
    
    return Solout;
    
}



