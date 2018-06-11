
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
#include "PostProcessTemplate.h"
#include "PostProcess.h"


void OneElementMatrix();
//Test: Buildconectivity and AutoBuild
void AutoBuildTest();

//ForceFunction
void FoFunc(const VecDouble &x, VecDouble &f);

//Test TwoElement Global Matrix and Analysis
void TestGlobalMatrix();
void solexact1(const VecDouble &co, VecDouble &result, Matrix &deriv);
void solexact2(const VecDouble &co, VecDouble &result, Matrix &deriv);
void solexact3(const VecDouble &co, VecDouble &result, Matrix &deriv);

GeoMesh *CreateGeoMesh(int nx, int ny, int dim, double lx, double ly);
int main ()
{
  //  OneElementMatrix();
 //   AutoBuildTest();
  //   TestGlobalMatrix();
    int n = 64;
   GeoMesh *gmesh = CreateGeoMesh(n,n,2,1,1);
   gmesh->BuildConnectivity();

    CompMesh *cmesh = new CompMesh(gmesh);
    cmesh->SetNumberMath(5);
     Matrix perm(2, 2, 1.);
    perm(1,0) = 0.;
    perm(0,1) = 0.;
    Poisson * pos = new Poisson(1, perm);
    pos->SetForceFunction(FoFunc);
    std::vector<MathStatement *> mathvec(gmesh->NumElements());


    Matrix projection(2,2,0.);
    projection(0,0)=1;
     projection(1,1)=1;
    MathStatement * bc0 = new L2Projection(-1, projection);
    bc0->SetExact(&solexact3);
    MathStatement * bc1 = new L2Projection(-2, projection);
     bc1->SetExact(&solexact1);
    MathStatement * bc2 = new L2Projection(-3, projection);
     bc2->SetExact(&solexact2);
    MathStatement * bc3 = new L2Projection(-4, projection);
     bc3->SetExact(&solexact1);
    cmesh->SetNumberMath(20);

    for(int iel=0; iel<gmesh->NumElements(); iel++){
        GeoElement *geoel = gmesh->Element(iel);
        if(geoel->Material()==1){
            mathvec[iel]=pos;

        }
        if(geoel->Material()==-1){
            mathvec[iel]=bc0;
            cmesh->SetMathStatement(iel, bc0);
        }
        if(geoel->Material()==-2){
            mathvec[iel]=bc1;

        }
        if(geoel->Material()==-3){
          mathvec[iel]=bc2;
        }
        if(geoel->Material()==-4){
            mathvec[iel]=bc3;
        }
    }

    cmesh->SetMathVec(mathvec);
    cmesh->SetDefaultOrder(1);
    cmesh->AutoBuild();
    for(int i=0; i<gmesh->NumElements(); i++){
        cmesh->GetElement(i)->SetStatement(mathvec[i]);
    }
   
    Analysis an(cmesh);

    an.RunSimulation();
    PostProcess *solpos = new PostProcessTemplate<Poisson>(&an);
    
    
    solpos->AppendVariable("Sol");
//      solpos->AppendVariable("DSol");
//        solpos->AppendVariable("Sol_exact");
//        solpos->AppendVariable("Force");
    
    an.PostProcessSolution("SolutionPost.vtk", *solpos);
    

    
//    VTKGeoMesh::PrintGMeshVTK(gmesh, "TESTVMALLA.vtk");
//    gmesh->Print(std::cout);
//    std::cout<<"bingo";
    
    
    return 0;
};
void solexact1(const VecDouble &co, VecDouble &result, Matrix &deriv)
{
    result[0] = 0;
    result[1] = 0;
    result[2] = 0.;
}
void solexact2(const VecDouble &co, VecDouble &result, Matrix &deriv)
{
    result[0] = 0;
    result[1] = 0;
    result[2] = 0.;
}
void solexact3(const VecDouble &co, VecDouble &result, Matrix &deriv)
{
    result[0] = 0;
    result[1] = 0;
    result[2] = 0.;
}
void FoFunc(const VecDouble &x, VecDouble &f){
    f.resize(2);
    
    double xP = x[0];
    double yP = x[1];
    
    
    double fx =  2*M_PI*M_PI*sin(M_PI*xP)*sin(M_PI*yP);
    double fy =  yP;
    
    f[0] = fx;
    f[1] = 0;
}
GeoMesh *CreateGeoMesh(int nx, int ny, int dim, double lx, double ly){
   
    double eps = 1.0e-10; // closed point tolerance
    double hx=lx/nx;
    double hy=ly/ny;
    int NumNodes =(nx+1)*(ny+1);
    int NumElements = nx*ny;
    Matrix pointsNiv(ny+1,nx+1,0);
    VecInt nodes;
    Matrix nodeCord(2,NumNodes,0);
    int npx = nx+1;
    int npy= ny+1;
    int n_node;
    GeoMesh *gmesh = new GeoMesh;
    gmesh->SetNumNodes(NumNodes);
    gmesh->SetNumElements(NumElements);
    VecDouble co(3,0);
    VecInt Nodes(4,0);
    //create nodes
    for(int j=0; j<npy; j++){
        for(int i=0; i<npx; i++){
            n_node = i + (j)*npx;
            nodes.push_back(n_node);
            pointsNiv(j,i)=n_node;
            nodeCord(0,i +(j)*npx) = hx * i;
            nodeCord(1,i +(j)*npx) = hy * j;
            co[0]=hx * i;
            co[1]=hy * j;
            gmesh->Node(n_node).SetCo(co);
        }
    }
    //create GeoElements elements
    int index;
    int MatId;
    for(int j=0; j<ny; j++){
        for(int i=0; i<nx; i++){
            Nodes[0]= i + (nx+1)*(j);
            Nodes[1]= i + 1 + (nx+1)*(j);
            Nodes[2]= i + 1 + (nx+1)*(j+1) ;
            Nodes[3]= i + (nx+1)*(j+1);
            index = i + j*nx;
            MatId=1;
            GeoElement *gel = new GeoElementTemplate<GeomQuad>(Nodes,MatId,gmesh,index);
            gmesh->SetElement(index, gel);
        }
    }
    //Boundary Conditions
    
    int bcL = -1;
    int bcT = -2;
    int bcR = -3;
    int bcB = -4;
    
    VecDouble x_min_max(2,0);
    VecDouble y_min_max(2,0);
    x_min_max[0] = 0.0;
    x_min_max[1] = lx;
    y_min_max[0] = 0.0;
    y_min_max[1] = ly;
    
    std::vector<std::pair<int, int > > set_bcL;
    std::vector<std::pair<int, int > > set_bcT;
    std::vector<std::pair<int, int > > set_bcR;
    std::vector<std::pair<int, int > > set_bcB;
    
    std::pair<int, int > chunk;
    
    double coord_val;
    for(int iel=0; iel<NumElements; iel++){
        GeoElement *gel = gmesh->Element(iel);
        int nnodes= gel->NNodes();
        
        for (int i_node = 0; i_node < nnodes ; i_node++) {
            
            int glob_n_index = gel->NodeIndex(i_node);
            coord_val = gmesh->Node(glob_n_index).Coord(0);
            
            if (fabs(x_min_max[0] - coord_val) < eps) {
                chunk.first = gel->GetIndex();
                chunk.second = glob_n_index;
                set_bcL.push_back(chunk);
            }
            
            if (fabs(x_min_max[1] - coord_val) < eps) {
                chunk.first = gel->GetIndex();
                chunk.second = glob_n_index;
                set_bcR.push_back(chunk);
            }

            if (dim == 2) {
                coord_val = gmesh->Node(gel->NodeIndex(i_node)).Coord(1);
                
                if (fabs(y_min_max[0] - coord_val) < eps) {
                    chunk.first = gel->GetIndex();
                    chunk.second = glob_n_index;
                    set_bcB.push_back(chunk);
                }
                
                if (fabs(y_min_max[1] - coord_val) < eps) {
                    chunk.first = gel->GetIndex();
                    chunk.second = glob_n_index;
                    set_bcT.push_back(chunk);
                }
            }
        }
    }
    
//    return gmesh;
    
    Nodes.resize(2);
    int max_gel_index;
    
    // Inserting Left elements
    for (int ibc = 0; ibc < set_bcL.size() - 1; ibc++) {
        if (set_bcL[ibc].first == set_bcL[ibc+1].first) {
            Nodes[0]= set_bcL[ibc].second;
            Nodes[1]= set_bcL[ibc+1].second;
            max_gel_index = gmesh->NumElements();
            gmesh->SetNumElements(max_gel_index+1);
            GeoElement *gel = new GeoElementTemplate<Geom1d>(Nodes,bcL,gmesh,max_gel_index);
            gmesh->SetElement(max_gel_index, gel);

        }
    }
    
    // Inserting Top elements
    for (int ibc = 0; ibc < set_bcT.size() - 1; ibc++) {
        if (set_bcT[ibc].first == set_bcT[ibc+1].first) {
            Nodes[0]= set_bcT[ibc].second;
            Nodes[1]= set_bcT[ibc+1].second;
            max_gel_index = gmesh->NumElements();
            gmesh->SetNumElements(max_gel_index+1);
            GeoElement *gel = new GeoElementTemplate<Geom1d>(Nodes,bcT,gmesh,max_gel_index);
            gmesh->SetElement(max_gel_index, gel);
            
        }
    }
    
    // Inserting Right elements
    for (int ibc = 0; ibc < set_bcR.size() - 1; ibc++) {
        if (set_bcR[ibc].first == set_bcR[ibc+1].first) {
            Nodes[0]= set_bcR[ibc].second;
            Nodes[1]= set_bcR[ibc+1].second;
            max_gel_index = gmesh->NumElements();
            gmesh->SetNumElements(max_gel_index+1);
            GeoElement *gel = new GeoElementTemplate<Geom1d>(Nodes,bcR,gmesh,max_gel_index);
            gmesh->SetElement(max_gel_index, gel);
            
        }
    }
    
    // Inserting Bottom elements
    for (int ibc = 0; ibc < set_bcB.size() - 1; ibc++) {
        if (set_bcB[ibc].first == set_bcB[ibc+1].first) {
            Nodes[0]= set_bcB[ibc].second;
            Nodes[1]= set_bcB[ibc+1].second;
            max_gel_index = gmesh->NumElements();
            gmesh->SetNumElements(max_gel_index+1);
            GeoElement *gel = new GeoElementTemplate<Geom1d>(Nodes,bcB,gmesh,max_gel_index);
            gmesh->SetElement(max_gel_index, gel);
            
        }
    }
    
    return gmesh;
    
}
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
    
 
    Analysis an(cmesh);
    an.RunSimulation();
   
    PostProcess *solpos = new PostProcessTemplate<Poisson>(&an);
   
    
    solpos->AppendVariable("Sol");
    //  solpos->AppendVariable("DSol");
//    solpos->AppendVariable("Sol_exact");
//    solpos->AppendVariable("Force");
    
    an.PostProcessSolution("SolutionPost.vtk", *solpos);
    
    
 
    
//    std::cout<<"\n"<<"Num Eq: ";
//    std::cout<<Assem.NEquations();
//    std::cout<<"\n";
//    Matrix globmat;
//    Matrix rhs;
//    Assem.Compute(globmat, rhs);
//    globmat.PrintM();
//    rhs.PrintM();
    
//    CompElement * cel = cmesh->GetElement(0);
//    cel->SetStatement(material);
//    Matrix ek;
//    Matrix ef;
//    cel->CalcStiff(ek, ef);
    
    
    
    
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
    
    cmesh.SetDefaultOrder(1);
    
    
    cmesh.AutoBuild();
    std::cout<<"\n"<<cmesh.GetNumberDOF()<<"\n";
    
    //    VTKGeoMesh::PrintGMeshVTK(&geomesh, "TESTVMALLA.vtk");
    //    geomesh.Print(std::cout);
    //    std::cout<<"bingo";
    //
    gmesh.Print(std::cout);

};

