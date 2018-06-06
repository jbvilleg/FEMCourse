
//

#include <iostream>
#include "GeomQuad.h"
#include "Geom1d.h"
#include "GeomTetrahedron.h"
#include "GeomTriangle.h"
#include "GeoElement.h"
#include "GeoElementTemplate.h"
#include "Poisson.h"
#include "Assemble.h"
#include "Topology1d.h"
#include "TopologyTriangle.h"
#include "TopologyQuad.h"
#include "TopologyTetrahedron.h"
#include "IntRule.h"
#include "IntRule1d.h"
#include "IntRuleQuad.h"
#include "IntRuleTetrahedron.h"
#include "IntRuleTriangle.h"
#include "DataTypes.h"
#include "VTKGeoMesh.h"
#include "ReadGmsh.h"
#include <math.h>
#include "Poisson.h"
#include "tpanic.h"
#include <functional>
#include "Analysis.h"
#include "Assemble.h"


void OneElementMatrix();
//Test: Buildconectivity and AutoBuild
void AutoBuildTest();

//ForceFunction
void FoFunc(const VecDouble &x, VecDouble &f);

//Test TwoElement Global Matrix
void TestGlobalMatrix();
int main ()
{
 //   OneElementMatrix();
 //   AutoBuildTest();
    TestGlobalMatrix();
    
    return 0;
};

void TestGlobalMatrix(){
    //Creating the geometric mesh
    GeoMesh *gmesh = new GeoMesh;

    gmesh->SetNumNodes(6);
    VecDouble co(3);
    co[0]=0;
    co[1]=0;
    co[2]=0;
    gmesh->Node(0).SetCo(co);
    co[0]=2;
    co[1]=0;
    co[2]=0;
    gmesh->Node(1).SetCo(co);
    co[0]=4;
    co[1]=0;
    co[2]=0;
    gmesh->Node(2).SetCo(co);
    co[0]=4;
    co[1]=2;
    co[2]=0;
    gmesh->Node(3).SetCo(co);
    co[0]=2;
    co[1]=2;
    co[2]=0;
    gmesh->Node(4).SetCo(co);
    co[0]=0;
    co[1]=2;
    co[2]=0;
    gmesh->Node(5).SetCo(co);

    
    Matrix perm(2,2,0.);
    perm(0,0)=1.;
    perm(1,1)=1.;
    
    int MatId1=1;
    int index1=0;
    int MatId2=2;
    int index2=1;
    VecInt  nodes1(4);
    nodes1[0]=0;
    nodes1[1]=1;
    nodes1[2]=4;
    nodes1[3]=5;
    VecInt  nodes2(4);
    nodes2[0]=1;
    nodes2[1]=2;
    nodes2[2]=3;
    nodes2[3]=4;
    
    GeoElement *gel1 = new GeoElementTemplate<GeomQuad>(nodes1,MatId1,gmesh,index1);
    GeoElement *gel2 = new GeoElementTemplate<GeomQuad>(nodes2,MatId2,gmesh,index2);
    
    gmesh->SetNumElements(2);
    gmesh->SetElement(index1, gel1);
    gmesh->SetElement(index2, gel2);
    gmesh->BuildConnectivity();
    gmesh->Print(std::cout);
    
    //Creating the computational Mesh
    CompMesh *cmesh = new CompMesh;
    Poisson *material = new Poisson(1,perm);
    material->SetForceFunction(FoFunc);
    std::vector<MathStatement *> mathvec(2);
    mathvec[0] = material;
    mathvec[1] = material;
    cmesh->SetMathVec(mathvec);
    
    cmesh->SetGeoMesh(gmesh);
    cmesh->SetNumberMath(2);
    cmesh->SetMathStatement(0, material);
    cmesh->SetMathStatement(1, material);
    cmesh->SetDefaultOrder(1);
        cmesh->AutoBuild();
    CompElement * cel1 = cmesh->GetElement(0);
    cel1->SetStatement(material);
    CompElement * cel2 = cmesh->GetElement(1);
    cel2->SetStatement(material);
    
    Assemble Assem(cmesh);
    std::cout<<"\n"<<"Num Eq: ";
    std::cout<<Assem.NEquations();
    Matrix globmat;
    Matrix rhs;
    Assem.Compute(globmat, rhs);
    
//    CompElement * cel = cmesh->GetElement(0);
//    cel->SetStatement(material);
//    Matrix ek;
//    Matrix ef;
//    cel->CalcStiff(ek, ef);
    
    
    
    
}


void FoFunc(const VecDouble &x, VecDouble &f){
    f.resize(2);
    
    double xP = x[0];
    double yP = x[1];
  
    
    double fx =  xP;
    double fy =  yP;
    
    f[0] = fx;
    f[1] = fy;
}


void OneElementMatrix(){
    //create a elemento geometrico.
    
    //Creating the geometric mesh
    GeoMesh *gmesh = new GeoMesh;
    int MatId=1;
    int index=0;
    VecInt  nodes(4);
    nodes[0]=0;
    nodes[1]=1;
    nodes[2]=2;
    nodes[3]=3;
    gmesh->SetNumNodes(4);
    VecDouble co(3);
    co[0]=0;
    co[1]=0;
    co[2]=0;
    gmesh->Node(0).SetCo(co);
    co[0]=2;
    co[1]=0;
    co[2]=0;
    gmesh->Node(1).SetCo(co);
    co[0]=2;
    co[1]=2;
    co[2]=0;
    gmesh->Node(2).SetCo(co);
    co[0]=0;
    co[1]=2;
    co[2]=0;
    gmesh->Node(3).SetCo(co);
    
    Matrix perm(2,2,0.);
    perm(0,0)=1.;
    perm(1,1)=1.;
    
    GeoElement *gel = new GeoElementTemplate<GeomQuad>(nodes,MatId,gmesh,index);
    gmesh->SetNumElements(1);
    gmesh->SetElement(index, gel);
    gmesh->BuildConnectivity();
    
    //Creating the computational Mesh
    CompMesh *cmesh = new CompMesh;
    Poisson *material = new Poisson(1,perm);
    material->SetForceFunction(FoFunc);
    std::vector<MathStatement *> mathvec(1);
    mathvec[0] = material;
    cmesh->SetMathVec(mathvec);
    
    cmesh->SetGeoMesh(gmesh);
    cmesh->SetNumberMath(1);
    cmesh->SetMathStatement(0, material);
    cmesh->SetDefaultOrder(1);
    cmesh->AutoBuild();
    
    CompElement * cel = cmesh->GetElement(0);
    cel->SetStatement(material);
    Matrix ek;
    Matrix ef;
    cel->CalcStiff(ek, ef);
    std::cout<<"\n";
    std::cout<<"MATRIX OF ELEMENT"<<"\n";
    
    ek.Print();
      std::cout<<"\n";
    std::cout<<"FORCE OF ELEMENT"<<"\n";
      std::cout<<"\n";
    ef.Print();
    
    std::cout<<"bingo"<<"\n";

    
};
void AutoBuildTest(){
    
    //Test: Buildconectivity and AutoBuild
    ReadGmsh read;
    GeoMesh gmesh;
    CompMesh cmesh;
    
    read.Read(gmesh, "TESTV.msh");
    cmesh.SetGeoMesh(&gmesh);
    Matrix perm(2,2,0.);
    perm(0,0)=1.;
    perm(1,1)=1.;
    
    int nel= gmesh.NumElements();
    for (int iel = 0; iel<nel; iel++) {
        int geoMatID = 1;
        
        // Materiais internos (Poisson)
        cmesh.SetNumberMath(iel+1);
        Poisson *material = new Poisson(geoMatID,perm);
        //  material->SetForceFunction(F_source);
        cmesh.SetMathStatement(iel, material);
        
    }
    
    cmesh.SetDefaultOrder(0);
    
    
    cmesh.AutoBuild();
    std::cout<<"\n"<<cmesh.GetNumberDOF()<<"\n";
    
    //    VTKGeoMesh::PrintGMeshVTK(&geomesh, "TESTVMALLA.vtk");
    //    geomesh.Print(std::cout);
    //    std::cout<<"bingo";
    //
    gmesh.Print(std::cout);

};

