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
   
    
    detjac = 0.0;
    int nrows = gradx.Rows();
    int ncols = gradx.Cols();
    int dim   = ncols;
    
    switch (dim) {
        case 0:
        {
            jac.Resize(dim,dim);
            axes.Resize(dim,3);
            jacinv.Resize(dim,dim);
            detjac=1.;
            
            break;
        }
            
        case 1:
        {
            jac.Resize(dim,dim);
            axes.Resize(dim,3);
            jacinv.Resize(dim,dim);
            jac.Zero();
            
            /**  Definitions: v1 -> is the xi_direction of the Gradient */
            VecDouble v_1(3,0.);
            
            for (int i = 0; i < nrows; i++) {
                v_1[i]  = gradx.GetVal(i,0);
            }
            
            double norm_v_1 = 0.;
            for(int i = 0; i < nrows; i++) {
                norm_v_1 += v_1[i]*v_1[i];
            }
            
            norm_v_1    = sqrt(norm_v_1);
            jac(0,0)    = norm_v_1;
            detjac      = norm_v_1;
            jacinv(0,0) = 1.0/detjac;
            
            detjac = fabs(detjac);
            
            for(int i=0; i < 3; i++) {
                axes(0,i) = v_1[i]/norm_v_1;
            }
            
        }
            break;
        case 2:
        {
            
            jac.Resize(dim,dim);
            axes.Resize(dim,3);
            jacinv.Resize(dim,dim);
            jac.Zero();
            
            /**  Definitions: v1 -> is the xi_direction of the Gradient, v2 -> is the eta_direction of the Gradient */
            VecDouble v_1(3,0.), v_2(3,0.);
            
            /**  Definitions: v_1_til and v_2_til -> asscoiated orthonormal vectors to v_1 and v_2 */
            VecDouble v_1_til(3,0.), v_2_til(3,0.);
            
            for (int i = 0; i < nrows; i++) {
                v_1[i]  = gradx.GetVal(i,0);
                v_2[i]  = gradx.GetVal(i,1);
            }
            
            double norm_v_1_til = 0.0;
            double norm_v_2_til = 0.0;
            double v_1_dot_v_2  = 0.0;
            
            for(int i = 0; i < 3; i++) {
                norm_v_1_til    += v_1[i]*v_1[i];
                v_1_dot_v_2     += v_1[i]*v_2[i];
            }
            norm_v_1_til = sqrt(norm_v_1_til);
            
            for(int i=0 ; i < 3; i++) {
                v_1_til[i]          = v_1[i] / norm_v_1_til; // Normalizing
                v_2_til[i]          = v_2[i] - v_1_dot_v_2 * v_1_til[i] / norm_v_1_til;
                norm_v_2_til   += v_2_til[i]*v_2_til[i];
            }
            norm_v_2_til = sqrt(norm_v_2_til);
            
            
            jac(0,0) = norm_v_1_til;
            jac(0,1) = v_1_dot_v_2/norm_v_1_til;
            jac(1,1) = norm_v_2_til;
            
            detjac = jac(0,0)*jac(1,1)-jac(1,0)*jac(0,1);
            
            jacinv(0,0) = +jac(1,1)/detjac;
            jacinv(1,1) = +jac(0,0)/detjac;
            jacinv(0,1) = -jac(0,1)/detjac;
            jacinv(1,0) = -jac(1,0)/detjac;
            
            detjac = fabs(detjac);
            
            
            for(int i=0; i < 3; i++) {
                v_2_til[i] /= norm_v_2_til; 
                axes(0,i)  = v_1_til[i];
                axes(1,i)  = v_2_til[i];
            }
            
        }
            break;
        case 3:
        {
            jac.Resize(dim,dim);
            axes.Resize(dim,3);
            jacinv.Resize(dim,dim);
            jac.Zero();
            
            for (int i = 0; i < nrows; i++) {
                jac(i,0)  = gradx.GetVal(i,0);
                jac(i,1)  = gradx.GetVal(i,1);
                jac(i,2)  = gradx.GetVal(i,2);
            }
            
            detjac -= jac(0,2)*jac(1,1)*jac(2,0);//- a02 a11 a20
            detjac += jac(0,1)*jac(1,2)*jac(2,0);//+ a01 a12 a20
            detjac += jac(0,2)*jac(1,0)*jac(2,1);//+ a02 a10 a21
            detjac -= jac(0,0)*jac(1,2)*jac(2,1);//- a00 a12 a21
            detjac -= jac(0,1)*jac(1,0)*jac(2,2);//- a01 a10 a22
            detjac += jac(0,0)*jac(1,1)*jac(2,2);//+ a00 a11 a22
            
            jacinv(0,0) = (-jac(1,2)*jac(2,1)+jac(1,1)*jac(2,2))/detjac;//-a12 a21 + a11 a22
            jacinv(0,1) = ( jac(0,2)*jac(2,1)-jac(0,1)*jac(2,2))/detjac;//a02 a21 - a01 a22
            jacinv(0,2) = (-jac(0,2)*jac(1,1)+jac(0,1)*jac(1,2))/detjac;//-a02 a11 + a01 a12
            jacinv(1,0) = ( jac(1,2)*jac(2,0)-jac(1,0)*jac(2,2))/detjac;//a12 a20 - a10 a22
            jacinv(1,1) = (-jac(0,2)*jac(2,0)+jac(0,0)*jac(2,2))/detjac;//-a02 a20 + a00 a22
            jacinv(1,2) = ( jac(0,2)*jac(1,0)-jac(0,0)*jac(1,2))/detjac;//a02 a10 - a00 a12
            jacinv(2,0) = (-jac(1,1)*jac(2,0)+jac(1,0)*jac(2,1))/detjac;//-a11 a20 + a10 a21
            jacinv(2,1) = ( jac(0,1)*jac(2,0)-jac(0,0)*jac(2,1))/detjac;//a01 a20 - a00 a21
            jacinv(2,2) = (-jac(0,1)*jac(1,0)+jac(0,0)*jac(1,1))/detjac;//-a01 a10 + a00 a11
            
            detjac = fabs(detjac);
            
            axes.Zero();
            axes(0,0) = 1.0;
            axes(1,1) = 1.0;
            axes(2,2) = 1.0;
            
        }
            break;
            
        default:
        {
            std::cout << " Object with wrong dimensions, unable to compute jacobian matrix. Dimension = " << dim << std::endl;
            DebugStop();
        }
            break;
    }
    

    
    
    
    
    
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
