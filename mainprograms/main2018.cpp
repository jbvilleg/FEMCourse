

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
#include "tpanic.h"

#include "Analysis.h"



int main ()
{

    ReadGmsh read;
    GeoMesh geomesh;
    // lectura de malla con cuatro elementos.!
    read.Read(geomesh, "TESTV.msh");
    
    VTKGeoMesh::PrintGMeshVTK(&geomesh, "TESTVMALLA.vtk");
//    geomesh.Print(std::cout);
//    std::cout<<"bingo";
//    
    geomesh.Print(std::cout);
    
    return 0;
};
