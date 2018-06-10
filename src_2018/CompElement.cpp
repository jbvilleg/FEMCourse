//
//  GeoMesh.cpp
//  FemSC
//
//  Created by Philippe Devloo on 16/04/18.
//

#include "CompElement.h"
#include "GeoElement.h"
#include "CompMesh.h"
#include "CompElementTemplate.h"
#include "DataTypes.h"
#include "MathStatement.h"
#include "PostProcessTemplate.h"
#include "PostProcess.h"

    // Default constructor of CompElement
CompElement::CompElement(){
    geoel = 0;
    *intrule = 0;
    mat = 0;
    *compmesh = 0;
    index = -1;
 
}
// Constructor of CompElement
CompElement::CompElement(int64_t ind, CompMesh *cmesh, GeoElement *geo){
    index = ind;
    compmesh = cmesh;
    geoel = geo;
    
}

// Copy constructor of CompElement
CompElement::CompElement(const CompElement &copy){
    index = copy.index;
    compmesh = copy.compmesh;
    geoel = copy.geoel;
    intrule = copy.intrule;
    mat = copy.mat;
    
}

// Operator of copy
CompElement &CompElement::operator=(const CompElement &copy){
   
    index = copy.index;
    compmesh = copy.compmesh;
    geoel = copy.geoel;
    intrule = copy.intrule;
    mat = copy.mat;
    return *this;

}

// Destructor of CompElement
CompElement::~CompElement(){
    
}

// Method for creating a copy of the element
CompElement *CompElement::Clone() const{
    
}

// Return the material object associated with the element
MathStatement *CompElement::GetStatement() const{
    return mat;
}

// Set the material object associated with the element
void CompElement::SetStatement(MathStatement *statement){
    mat=statement;
}

// Return integration rule established
IntRule *CompElement::GetIntRule() const{
    return intrule;
}

// Set integration rule established
void CompElement::SetIntRule(IntRule *intrul){
    intrule =intrul;
}


// Set element index
void CompElement::SetIndex(int64_t ind){
    index=ind;
}

// Return the geometric element associated
GeoElement *CompElement::GetGeoElement() const{
    return geoel;
}

// Set the geometric element associated
void CompElement::SetGeoElement(GeoElement *element){
    geoel = element;
}

// Return a pointer to the element computational mesh
CompMesh *CompElement::GetCompMesh() const{
    return compmesh;
    
}

// Set a pointer to the element computational mesh
void CompElement::SetCompMesh(CompMesh *mesh){
    compmesh = mesh;
}

// Initialize integration points data object
void CompElement::InitializeIntPointData(IntPointData &data) const{
    int dim = this->Dimension();
    int nshape = this->NShapeFunctions();
    int nstate = this->GetStatement()->NState();
    data.phi.resize(nshape,1);
    data.dphidksi.Resize(dim,nshape);
    data.dphidx.Resize(dim,nshape);
    data.axes.Resize(dim,3);
  
    data.x.resize(3);
    data.ksi.resize(dim);
    data.solution.resize(nstate);
    data.dsoldksi.Resize(dim,nstate);
    data.dsoldx.Resize(dim,nstate);
 
    
}

// Compute and fill integration points data object
void CompElement::ComputeRequiredData(IntPointData &data, VecDouble &intpoint) const{
    
    GeoElement *gel = this->GetGeoElement();
    Matrix gradx,Jac,JacInv;
    gel->X(intpoint, data.x);
    gel->GradX(intpoint, data.x, data.gradx);
    gel->Jacobian(data.gradx, Jac, data.axes, data.detjac, JacInv);
    this->ShapeFunctions(intpoint, data.phi, data.dphidksi);
    this->Convert2Axes(data.dphidksi, JacInv, data.dphidx);

}
void CompElement::Convert2Axes(const Matrix &dphi, const Matrix &jacinv, Matrix &dphidx) const{
    int nshape = dphi.Cols();
    int dim = dphi.Rows();
    dphidx.Resize(dim,nshape);
    int ieq;
    switch(dim){
        case 0:
        {
            
        }
            break;
        case 1:
        {
            dphidx = dphi;
            dphidx = dphidx*jacinv.GetVal(0,0);
        }
            break;
        case 2:
        {
            for(ieq = 0; ieq < nshape; ieq++) {
                dphidx(0,ieq) = jacinv.GetVal(0,0)*dphi.GetVal(0,ieq) + jacinv.GetVal(1,0)*dphi.GetVal(1,ieq);
                dphidx(1,ieq) = jacinv.GetVal(0,1)*dphi.GetVal(0,ieq) + jacinv.GetVal(1,1)*dphi.GetVal(1,ieq);
            }
        }
            break;
        case 3:
        {
            for(ieq = 0; ieq < nshape; ieq++) {
                dphidx(0,ieq) = jacinv.GetVal(0,0)*dphi.GetVal(0,ieq) + jacinv.GetVal(1,0)*dphi.GetVal(1,ieq) + jacinv.GetVal(2,0)*dphi.GetVal(2,ieq);
                dphidx(1,ieq) = jacinv.GetVal(0,1)*dphi.GetVal(0,ieq) + jacinv.GetVal(1,1)*dphi.GetVal(1,ieq) + jacinv.GetVal(2,1)*dphi.GetVal(2,ieq);
                dphidx(2,ieq) = jacinv.GetVal(0,2)*dphi.GetVal(0,ieq) + jacinv.GetVal(1,2)*dphi.GetVal(1,ieq) + jacinv.GetVal(2,2)*dphi.GetVal(2,ieq);
            }
        }
            break;
        default:
        {
            std::cout << "Error at " << __PRETTY_FUNCTION__ << " please implement the " << dim << "d Jacobian and inverse\n" <<std::endl;
        }
    }
    
}


// Compute the element stifness matrix and force vector
void CompElement::CalcStiff(Matrix &ek, Matrix &ef) const{
    
    MathStatement *material = GetStatement();
    if(!material){
        std::cout << " Material == NULL " << std::endl;
        return;
    }
    //this->InitializeElementMatrix(ek,ef);
    IntPointData data;
    this->InitializeIntPointData(data);
    double weight =1.;
    int nintrulepoints = GetIntRule()->NPoints();
    int dim = Dimension();
    VecDouble intpoint(dim,0.);
    //resize! ok
    int nshape = data.phi.size();
    int nstate = GetStatement()->NState();
    ek.Resize(nshape*nstate,nshape*nstate);
    ef.Resize(nshape*nstate,nstate);
    ek.Zero();
    ef.Zero();
    
    for (int intd_id = 0; intd_id < nintrulepoints; intd_id++) {
        intrule->Point(intd_id, data.ksi, data.weight);
        this->ComputeRequiredData(data, data.ksi);
        weight=data.weight;
        weight *=fabs(data.detjac);
        material->Contribute(data, weight, ek, ef);
    }

}
void CompElement::Solution(VecDouble &intpoint, int var, VecDouble &sol, TMatrix &dsol) const{
    
    IntPointData data;
    this->InitializeIntPointData(data);
    GetMultiplyingCoeficients(data.coefs);
    MathStatement *mat = this->GetStatement();
    
    ComputeRequiredData(data,intpoint);
    data.ComputeSolution();
    sol.resize(2);
    dsol.Resize(data.dsoldx.Rows(),data.dsoldx.Cols());
    
    sol=mat->PostProcessSolution(data,var);
    dsol=data.dsoldx;
    
}

