    //
//  Geom1d.h
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.

#include "ShapeQuad.h"
#include "tpanic.h"
#include "TopologyQuad.h"
#include "Topology1d.h"
#include "TopologyTriangle.h"

 int TopologyTriangle::NSideNodes(int side){
     
     if(side<3){
         return 1;
     }
     
     if (side>=3 && side <=5){
         return 2;
     }
     
     if (side==6){
         return 3;
     }
    
     if (side > 6){
         std::cout<<"El triangulo tiene 7 lados (de 0 a 6)";
         DebugStop();
     }
     
}

/// local node index of a node associated with a side
 int TopologyTriangle::SideNodeIndex(int side, int node){
     if (side <3 && node==0){
         return side;
     }
     
     if (side>=3 && side<6 && node <2){
         return (side-node)%3;
     }
     
     if (side==6 && node<3){
         return node;
     }
     
     
     if (side > 6){
         std::cout<<"El triangulo tiene 7 lados (0 a 6)";
         DebugStop();
     }
    
}

/// return the enumerated element type
 ElementType TopologyTriangle::Type(){
     return ElementType(2);
}
