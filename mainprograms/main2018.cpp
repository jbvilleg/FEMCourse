

#include <iostream>
#include <math.h>
#include "IntRule.h"
#include "IntRule1d.h"
#include "IntRuleQuad.h"
#include "IntRuleTetrahedron.h"
#include "IntRuleTriangle.h"
#include "Topology1d.h"
#include "TopologyTriangle.h"
#include "TopologyQuad.h"
#include "TopologyTetrahedron.h"
#include "Shape1d.h"
#include "DataTypes.h"
#include "GeoMesh.h"
#include "ReadGmsh.h"
#include "VTKGeoMesh.h"
#include "CompElementTemplate.h"
#include "GeoElement.h"
#include "CompElement.h"
#include "CompMesh.h"
#include "GeomQuad.h"
#include "GeoMesh.h"
#include "GeoElementTemplate.h"
#include "Poisson.h"
#include "CompElementTemplate.h"
#include "Shape1d.h"
#include "Geom1d.h"

using std::cout;
using std::endl;
using std::cin;

int main ()
{
    
    
    
    //-------------------------------------
    VecDouble vec1;
    CompMesh *compmesh;
    ReadGmsh read;
    
    GeoMesh geomesh;
    
    read.Read(geomesh, "TESTV.msh");
    
    
    CompElementTemplate<ShapeQuad> jos;
    Poisson *poasson;
    jos.SetGeoElement(geomesh.Element(12));
    jos.SetCompMesh(compmesh);
    IntRuleQuad *intrule = new IntRuleQuad();
    int order = 3;
    intrule->SetOrder(order);
    jos.SetIntRule(intrule);
    Matrix perm(2, 2);
    perm(1,0) = 0.;
    perm(0,1) = 0.;
    Matrix perm3;
    perm3=perm;
    Poisson *poasson2 = new Poisson(perm);
    jos.SetStatement(poasson2);
    
    Matrix ek;
    Matrix ef;
//    jos.CalcStiff(ek, ef);
//    compmesh.SetGeoMesh(&geomesh);
//    compmesh.AutoBuild();
    
    
    geomesh.BuildConnectivity();
    geomesh.Element(2);
    VTKGeoMesh::PrintGMeshVTK(&geomesh, "TESTVILLEGAS.vtk");
    geomesh.Print(cout);
    
    
    
    
    return 0;
    
}
