//
//  GeoMesh.cpp
//  FemSC
//
//  Created by Philippe Devloo on 16/04/18.
//

#include "GeoNode.h"
#include "GeoElement.h"
#include <string>
#include "GeoMesh.h"
#include "tpanic.h"
#include "GeoElementSide.h"


// Function to print results
void GeoNode::Print(std::ostream &out){
    out << "Node : fId = " << "?";
    out << "    Coordinates";
    for(int i=0;i<3;i++) out << "\t" << Coord(i);
    out << "\n";
}

