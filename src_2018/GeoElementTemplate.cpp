#include "GeoMesh.h"
#include "GeoElement.h"
#include "GeoElementTemplate.h"
#include "GeomTriangle.h"
#include "GeomTetrahedron.h"
#include "GeomQuad.h"
#include "Geom1d.h"
#include "tpanic.h"


template<class TGeom>
GeoElementTemplate<TGeom>::GeoElementTemplate(const VecInt &nodeindices, int materialid, GeoMesh *gmesh, int index) : GeoElement(materialid,gmesh,index){
    
    Geom.SetNodes(nodeindices);
    for(int iside=0;iside<Geom.nSides;iside++){
        Geom.SetNeighbour(iside,GeoElementSide(this,iside));
    }
    
}


template<class TGeom>
GeoElementTemplate<TGeom>::GeoElementTemplate(const GeoElementTemplate &copy) : GeoElement(copy){
    
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
void GeoElementTemplate<TGeom>::X(const VecDouble &xi, VecDouble &x){
    
    Matrix NodeCo(3,NNodes(),0.);
    for (int in=0; in<NNodes(); in++) {
        int nodid = this->NodeIndex(in);
        NodeCo(0,in)=GMesh->Node(nodid).Co()[0];
        NodeCo(1,in)=GMesh->Node(nodid).Co()[1];
        NodeCo(2,in)=GMesh->Node(nodid).Co()[2];
    }
    
    Geom.X(xi,NodeCo,x);

}

template<class TGeom>
void GeoElementTemplate<TGeom>::GradX(const VecDouble &xi, VecDouble &x, Matrix &gradx){

    Matrix NodeCo(3,NNodes(),0.);
    for (int in=0; in<NNodes(); in++) {
        int nodid = NodeIndex(in);
        NodeCo(0,in)=GMesh->Node(nodid).Co()[0];
        NodeCo(1,in)=GMesh->Node(nodid).Co()[1];
        NodeCo(2,in)=GMesh->Node(nodid).Co()[2];
    }
    Geom.GradX(xi, NodeCo, x, gradx);
    
    
}

template<class TGeom>
int GeoElementTemplate<TGeom>::WhichSide(VecInt &SideNodeIds){
    int nSidesEl = NSides();
    int nNodesSide = SideNodeIds.size();
    for(int side=0; side<nSidesEl; side++){
        
        //NodeIndex1D
        if(NSideNodes(side)==2 && nNodesSide==2){
            int nIndex1 = NodeIndex(SideNodeIndex(side, 0));
            int nIndex2 = NodeIndex(SideNodeIndex(side, 1));
            if((nIndex1 == SideNodeIds[0] && nIndex2 == SideNodeIds[1]) ||
               (nIndex2 == SideNodeIds[0] && nIndex1 == SideNodeIds[1])){
                return side;
            }
        }
        //revisar
        
//        if(NSideNodes(side)==3 && nNodesSide==3){
//            int nIndex1 = NodeIndex(SideNodeIndex(side, 0));
//            int nIndex2 = NodeIndex(SideNodeIndex(side, 1));
//            int nIndex3 = NodeIndex(SideNodeIndex(side, 2));
//            if((nIndex1 == SideNodeIds[0] && nIndex2 == SideNodeIds[1]) ||
//               (nIndex2 == SideNodeIds[0] && nIndex1 == SideNodeIds[1])){
//                return side;
//            }
//        }

    }

}


template<class TGeom>
void GeoElementTemplate<TGeom>::Jacobian(const Matrix &gradx, Matrix &jac,Matrix &axes, double &detjac, Matrix &jacinv){
   
    
    
    
}

template<class TGeom>
void GeoElementTemplate<TGeom>::Print(std::ostream &out){
    
    out << "ElType " << Type() << " matid " << MaterialId << " index " << GetIndex() << " nodes ";
    int i;
    for(i=0; i<NNodes(); i++) out << NodeIndex(i) << ' ';
    out << std::endl;
    
    for (i = 0;i < NSides();i++) {
        out << "Neighbours for side   " << i << " : ";
        GeoElementSide neighbour = Neighbour(i);
        GeoElementSide thisside(this,i);
        if(!(neighbour.Element()!=0&&neighbour.Side()>-1))
        {
            out << "No neighbour\n";
        }
        else {
            while ((neighbour == thisside)==false ){
                out << neighbour.Element()->GetIndex() << "/" << neighbour.Side() << ' ';
                neighbour = neighbour.Neighbour();
            }
            out << std::endl;
        }
        
    }
    
}

template class GeoElementTemplate<GeomTriangle>;
template class GeoElementTemplate<Geom1d>;
template class GeoElementTemplate<GeomQuad>;
template class GeoElementTemplate<GeomTetrahedron>;
