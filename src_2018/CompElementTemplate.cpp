//
//  CompElementTemplate.h
//  FemCourse
//
//  Created by Philippe Devloo on 24/04/18.
//

#include "CompElementTemplate.h"
#include "CompElement.h"
#include "GeoElementTemplate.h"
#include "GeoElement.h"
#include "Shape1d.h"
#include "ShapeQuad.h"
#include "ShapeTetrahedron.h"
#include "ShapeTriangle.h"
#include "CompMesh.h"
#include "tpanic.h"
#include "DataTypes.h"
#include "MathStatement.h"
#include "DOF.h"


template<class Shape>
CompElementTemplate<Shape>::CompElementTemplate() : CompElement(){
    dofindexes.resize(0);
    intrule=0;
}

template<class Shape>
CompElementTemplate<Shape>::CompElementTemplate(int64_t ind, CompMesh *cmesh,  GeoElement *geo) : CompElement(ind,cmesh,geo){
    

     int order = cmesh->GetDefaultOrder();
     intrule.SetOrder(order*2);
     SetIntRule(&intrule);

}

template<class Shape>
CompElementTemplate<Shape>::CompElementTemplate(const CompElementTemplate &copy) : CompElement(copy){
    dofindexes=copy.dofindexes;
    intrule=copy.intrule;
}

template<class Shape>
CompElementTemplate<Shape> &CompElementTemplate<Shape>::operator=(const CompElementTemplate &copy){
    dofindexes=copy.dofindexes;
    intrule=copy.intrule;
    return *this;
}

template<class Shape>
CompElementTemplate<Shape>::~CompElementTemplate(){
    
}

template<class Shape>
void CompElementTemplate<Shape>::SetNDOF(int64_t ndof){
    dofindexes.resize(ndof);
}

template<class Shape>
void CompElementTemplate<Shape>::SetDOFIndex(int i, int64_t dofindex){
    dofindexes[i]=dofindex;
}

template<class Shape>
int64_t CompElementTemplate<Shape>::GetDOFIndex(int i){
    return dofindexes[i];
}

template<class Shape>
CompElement * CompElementTemplate<Shape>::Clone() const{
    
}

template<class Shape>
void  CompElementTemplate<Shape>::ShapeFunctions(const VecDouble &intpoint, VecDouble &phi, Matrix &dphi) const{
    
    VecInt orders(NDOF());
    int ndof = NDOF();
    CompMesh *cmesh =this->GetCompMesh();
    for (int ic=0; ic<ndof; ic++) {
        orders[ic]=cmesh->GetDOF(ic).GetOrder();
    }
    Shape::Shape(intpoint, orders, phi, dphi);
    
}
template<class Shape>
void CompElementTemplate<Shape>::GetMultiplyingCoeficients(VecDouble &coefs) const{
    
    CompMesh *cmesh = this->GetCompMesh();
    int neq = cmesh->Solution().size();
    VecInt iGlob;
    coefs.resize(0);
    
    int ndofel = NDOF();
    int indexdof = 0;
    for (int idof =0; idof<ndofel; idof++) {
        int idcon = dofindexes[idof];
        DOF dof = cmesh->GetDOF(idcon);
        int nshape = dof.GetNShape();
        int nstat = dof.GetNState();
        for(int i=0; i<nshape*nstat; i++) {
            iGlob.resize(indexdof+1);
            coefs.resize(indexdof+1);
            iGlob[indexdof] = dof.GetFirstEquation()+i;
            coefs[indexdof] = cmesh->Solution()[iGlob[indexdof]];
            indexdof++;
        }
    }
    
}

template<class Shape>
int CompElementTemplate<Shape>::NShapeFunctions() const{
 
    int ndof = NDOF();
    int nshape = 0;
    for (int ic=0; ic<ndof; ic++) {
        nshape +=NShapeFunctions(ic);
    }
    return nshape;
}

template<class Shape>
int CompElementTemplate<Shape>::NDOF() const{
    return dofindexes.size();
}

/// returns the number of shape functions stored in the DOF data structure
template<class Shape>
int CompElementTemplate<Shape>::NShapeFunctions(int doflocindex) const{
    CompMesh cmesh = *this->GetCompMesh();
    DOF dofex = cmesh.GetDOF(doflocindex);
    return cmesh.GetDOF(doflocindex).GetNShape();
}

/// uses the Shape template class to compute the number of shape functions
template<class Shape>
int CompElementTemplate<Shape>::ComputeNShapeFunctions(int doflocindex, int order){
//    dofindexes.resize(doflocindex+1);
//    dofindexes[doflocindex]=doflocindex;
    return Shape::NShapeFunctions(doflocindex,order);
}




template class CompElementTemplate<Shape1d>;
template class CompElementTemplate<ShapeQuad>;
template class CompElementTemplate<ShapeTriangle>;
template class CompElementTemplate<ShapeTetrahedron>;
