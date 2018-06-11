//
//  L2Projection.h
//  FemCourse
//
//  Created by Philippe Devloo on 24/04/18.
//

#include "L2Projection.h"
#include "tpanic.h"

L2Projection::L2Projection(){
    
}

L2Projection::L2Projection(int materialid, Matrix &perm){
    SetMatID(materialid);
    projection=perm;
}

L2Projection::L2Projection(const L2Projection &copy){
    projection=copy.projection;
    forceFunction=copy.forceFunction;
}

L2Projection &L2Projection::operator=(const L2Projection &copy){
    projection=copy.projection;
    forceFunction=copy.forceFunction;
    return *this;
}

L2Projection *L2Projection::Clone() const{
    
}

L2Projection::~L2Projection(){
    
}

Matrix L2Projection::GetProjectionMatrix() const{
    return projection;
}

void L2Projection::SetProjectionMatrix(const Matrix &proj){
    projection=proj;
}

int L2Projection::NEvalErrors() const{
    return 3;
}

int L2Projection::NState() const{
    return 1;
}

int L2Projection::VariableIndex(const PostProcVar var) const{
    return 0;
}

// Return the variable index associated with the name
L2Projection::PostProcVar L2Projection::VariableIndex(const std::string &name){
    return ENone;
}

// Return the number of variables associated with the variable indexed by var. Param var Index variable into the solution, is obtained by calling VariableIndex
int L2Projection::NSolutionVariables(const PostProcVar var){
    return 0;
}

void L2Projection::Contribute(IntPointData &data, double weight, Matrix &EK, Matrix &EF) const{
    double BigNumber =MathStatement::gBigNumber;
    int nvars = this->NState();
    int nshape = data.phi.size();
    
    VecDouble result(data.x.size());
    Matrix deriv(data.x.size(), data.x.size());
    
     std::function<void (const VecDouble &loc, VecDouble &result, Matrix &deriv)> fExact = GetExact();
    if(fExact){
        fExact(data.x, result, deriv);
    }else{
        std::cout<<"revisar condiciones de contorno";
        DebugStop();
    }

    
    for (int i = 0; i < nshape; i++) {
        EF(i, 0) += weight * data.phi[i] * result[0] * BigNumber;
        for (int j = 0; j < nshape; j++) {
            EK(i, j) += weight * data.phi[i] * data.phi[j]*BigNumber;
         }
    }
    
}

// Method to implement error over element's volume
void L2Projection::ContributeError(IntPointData &integrationpointdata, VecDouble &u_exact, Matrix &du_exact, VecDouble &errors) const{
    
    return;
}


std::vector<double> L2Projection::PostProcessSolution(const IntPointData &integrationpointdata, const int var) const{
    VecDouble axl2(2,0.);
    return axl2;
}

