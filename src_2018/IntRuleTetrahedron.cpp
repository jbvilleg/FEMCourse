//
//  Geom1d.h
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.

#include "IntRuleQuad.h"
#include "IntRule1d.h"
#include "IntRuleTriangle.h"
#include "IntRuleTetrahedron.h"

IntRuleTetrahedron::IntRuleTetrahedron(){
    fOrder=0;
    SetOrder(fOrder);
}

IntRuleTetrahedron::IntRuleTetrahedron(int order){
    fOrder=order;
    SetOrder(fOrder);
}

void IntRuleTetrahedron::SetOrder(int order){
    fOrder=order;
    
    
    
    
    
    
    
    
}

