//
//  Geom1d.h
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//



#include "DataTypes.h"
#include "GeoElementSide.h"
#include "GeoElement.h"
#include "CompElement.h"
#include "CompMesh.h"
#include "CompElementTemplate.h"
#include "GeoElementTemplate.h"
#include "Shape1d.h"
#include "tpanic.h"


GeoElement::GeoElement(){
   
}


GeoElement::GeoElement(const GeoElement &copy){
    this->MaterialId=copy.MaterialId;
    this->GMesh = copy.GMesh;
}


CompElement *GeoElement::CreateCompEl(CompMesh *mesh, int64_t index){
    
    GeoElement *gel = mesh->GetGeoMesh()->Element(index);
    
    switch (gel->Type()) {
        case EOned:
            return new  CompElementTemplate<Shape1d>(index,mesh,gel);
            break;
        case EQuadrilateral:
            return new CompElementTemplate<ShapeQuad>(index,mesh,gel);
            break;
        case ETriangle:
            return new CompElementTemplate<ShapeTriangle>(index,mesh,gel);
            break;
        case EPiramide:
            DebugStop();
            break;
        case EPrisma:
            DebugStop();
            break;
        case ETetraedro:
//            return new CompElementTemplate<ShapeTetrahedron>(index,mesh,gel);
            break;
        default:
            DebugStop();
            break;
    }
    return 0;
    
}

 GeoElement::~GeoElement(){
     
}


void GeoElement::Print(std::ostream &out){
    
}

