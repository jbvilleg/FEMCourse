    //
//  Geom1d.h
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.

#include "ShapeQuad.h"
#include "tpanic.h"
#include "TopologyQuad.h"
#include "Topology1d.h"

/// Number of nodes associated with a side
int TopologyQuad::NSideNodes(int side){
   
    if(side>=0 && side<=3){
        return 1;
    }
    
    if(side>3 && side <8){
        return 2;
    }
    
    if(side==8){
        return 2;
    }
    if(side >= 9 ){
        //El Cuadrilatero tiene 9 lados (de 0 a 8)
        DebugStop();
    }
    
}

/// local node index of a node associated with a side
int TopologyQuad::SideNodeIndex(int side, int node){
    
    if((side>=0 && side<=3) && node==0){
        return side;
    }
    if((side>3 && side<=7) && node<2){
        return (side+node)%4;
    }
    if((side==8) && node<4){
        return node;
    }
    //Robustes del codigo VERIFICAR
    if (side>8 ){
        std::cout<<"El Cuadrilatero tiene 9 lados (de 0 a 8)<<"<<"\n";
        DebugStop();
    }
    if (node>4 ){
        std::cout<<"El Cuadrilatero tiene 4 nodos (de 0 a 3)<<"<<"\n";
        DebugStop();
    }
    
}

/// return the enumerated element type
ElementType TopologyQuad::Type(){
    
    return ElementType(3);
}
