//
//  Assemble.h
//  FemCourse
//
//  Created by Philippe Devloo on 08/05/18.

#include "DataTypes.h"
#include "CompMesh.h"
#include "PostProcess.h"
#include "Assemble.h"
#include "CompMesh.h"
#include "CompElement.h"
#include "Analysis.h"

Analysis::Analysis(): cmesh(0), Solution(), GlobalSystem(), RightHandSide(){
    
    
}

Analysis::Analysis(const Analysis &cp){
    cmesh=cp.cmesh;
    Solution=cp.Solution;
    GlobalSystem=cp.GlobalSystem;
    RightHandSide=cp.RightHandSide;
}

Analysis &Analysis::operator=(const Analysis &cp){
    cmesh=cp.cmesh;
    Solution=cp.Solution;
    GlobalSystem=cp.GlobalSystem;
    RightHandSide=cp.RightHandSide;
    return *this;
    
}

Analysis::Analysis(CompMesh *mesh) : cmesh(mesh){
    
}

void Analysis::SetMesh(CompMesh *mesh){
    cmesh=mesh;
}

CompMesh *Analysis::Mesh() const{
    return cmesh;
}

void Analysis::RunSimulation(){
    
    Assemble Assem(cmesh);
    int neq = Assem.NEquations();
    Matrix K(neq,neq);
    Matrix F(neq,1);

    
    Assem.Compute(K, F);
    K.Print();
    GlobalSystem = K;
    RightHandSide = F;
    K.Solve_LU(F);
    Solution=F;
    
    std::vector<double> lsol(Solution.Rows(),0.);
    for (int is=0; is<Solution.Rows(); is++) {
        lsol[is]=Solution(is,0);
    }
    cmesh->LoadSolution(lsol);

}

void Analysis::PostProcessSolution(const std::string &filename, PostProcess &defPostProc) const{
    
    
}
