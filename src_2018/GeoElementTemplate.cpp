

#include <stdio.h>
#include "GeoElementTemplate.h"
#include "GeoNode.h"
#include "Geom1d.h"
#include "GeomQuad.h"
#include "GeomTriangle.h"
#include "GeomTetrahedron.h"
#include "GeoElement.h"
#include "GeoMesh.h"
#include "tpanic.h"

/// constructor
template<class TGeom>
GeoElementTemplate<TGeom>::GeoElementTemplate(const VecInt &nodeindices, int materialid, GeoMesh *gmesh) : GeoElement(materialid,gmesh){
    Geom.SetNodes(nodeindices);
}


template<class TGeom>
GeoElementTemplate<TGeom>::GeoElementTemplate(const GeoElementTemplate &copy){
    Geom = copy.Geom;
}

template<class TGeom>
GeoElementTemplate<TGeom> &GeoElementTemplate<TGeom>::operator=(const GeoElementTemplate &copy){
    Geom = copy.Geom;
    return *this;
}

/// return the enumerated element type
template<class TGeom>
ElementType GeoElementTemplate<TGeom>::Type(){
    return TGeom::Type();
    
}

template<class TGeom>
void GeoElementTemplate<TGeom>::X(const VecDouble &xi,VecDouble &x){
    
    Matrix NodeCoor(3,NNodes());
    GeoNode np;
    
    int i,j;
    for(i=0; i<this->NNodes(); i++) {
        int id = this->NodeIndex(i);
        np = GMesh->Node(id);
        for(j=0;j<3;j++) {
            NodeCoor(j,i) = np.Coord(j);
        }
    }
    
    Geom.X(xi, NodeCoor, x);
    
    
}

template<class TGeom>
void GeoElementTemplate<TGeom>::GradX(const VecDouble &xi,VecDouble &x,Matrix &gradx){

    Matrix NdCor(3,this->NNodes(),0.);
    GeoNode np;
    
    int i,j;
    for(i=0; i<this->NNodes(); i++) {
        np = GMesh->Node(i);
        for(j=0;j<3;j++) {
            NdCor(j,i) = np.Coord(j);
        }
    }
    
    Geom.GradX(xi, NdCor, x, gradx);
    
}

template<class TGeom>
void GeoElementTemplate<TGeom>::Print(std::ostream &out){
    
    int i;
    out << "\n"<< "Nod Index: ";

    for(i=0; i<NNodes(); i++) {
        out << NodeIndex(i) << ' ';
    }
    
      out << std::endl;
    
    for (i = 0; i < NSides();i++)
    {
        out << "Neighbours for side   " << i << " : ";
        GeoElementSide neighbour = Neighbour(i);
        GeoElementSide thisside(this,i);
            if(!(neighbour.Element()!= 0 && neighbour.Side()>-1))
            {
                out << "No neighbour\n";
            }
     
    }
        
}

template class GeoElementTemplate<GeomTriangle>;
template class GeoElementTemplate<Geom1d>;
template class GeoElementTemplate<GeomQuad>;
//template class GeoElementTemplate<GeomTetrahedron>;
