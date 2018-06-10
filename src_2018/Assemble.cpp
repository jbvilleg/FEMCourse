//
//  Assemble.h
//  FemCourse
//
//  Created by Philippe Devloo on 08/05/18.
//


#include "DataTypes.h"
#include "CompMesh.h"
#include "PostProcess.h"
#include "Assemble.h"
#include "CompMesh.h"
#include "CompElement.h"
#include "GeoElement.h"
#include "GeoElementTemplate.h"

    
Assemble::Assemble(){
    
}
    
Assemble::Assemble(CompMesh *mesh){
    cmesh=mesh;
}
    
Assemble::Assemble(const Assemble &copy){
    cmesh=copy.cmesh;
}
    
Assemble &Assemble::operator=(const Assemble &copy){
    cmesh = copy.cmesh;
    return *this;
}
    
void Assemble::SetMesh(CompMesh *mesh){
    cmesh= mesh;
}
    
    /// Compute the total number of equations
int64_t Assemble::NEquations(){
    
    int64_t numEq = 0;
    int64_t numDoF = cmesh->GetNumberDOF();
    
    for (int idof=0; idof<numDoF; idof++) {
        DOF dof = cmesh->GetDOF(idof);
        int dofsize = dof.GetNShape()*dof.GetNState();
        numEq += dofsize;
    }
    return numEq;
    
}
    
    /// Optimize the bandwidth of the global system of equations
void Assemble::OptimizeBandwidth(){
    
}
    
    /// Compute the global stiffness matrix and right hand side
void Assemble::Compute(Matrix &globmat, Matrix &rhs){
 
    int neq = NEquations();
    int nel= cmesh->GetElementVec().size();
    globmat.Resize(neq, neq);
    globmat.Zero();
    globmat.Print();
    rhs.Resize(neq, 1);
    
    for (int el=0; el<nel; el++) {
        
        CompElement *cel=cmesh->GetElement(el);
        Matrix EK,EF;
        cel->CalcStiff(EK, EF);
        EK.Print();
        VecInt IndexG(neq,0);
        int ndofElm = cel->NDOF();
        int indexdof = 0;
        for (int idof =0; idof<ndofElm; idof++) {
            int idcon = cel->GetDOFIndex(idof);
            DOF dof = cmesh->GetDOF(idcon);
            int nshape = dof.GetNShape();
            int nstat = dof.GetNState();
            for(int i=0; i<nshape*nstat; i++) {
                IndexG[indexdof] = dof.GetFirstEquation()+i;
                indexdof++;
            }
        }
        
        for (int i=0; i<EK.Rows(); i++) {
            rhs(IndexG[i],0)+=EF(i,0);
            for (int j=0; j<EK.Rows(); j++) {
                globmat(IndexG[i],IndexG[j])+=EK(i,j);
            }
        }
//
    //    globmat.PrintM();
    }
}
