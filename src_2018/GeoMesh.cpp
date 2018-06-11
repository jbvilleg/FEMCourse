//  GeoMesh.cpp
//  FemSC
//
//  Created by Philippe Devloo on 16/04/18.


#include "GeoNode.h"
#include "GeoElement.h"
#include <string>
#include "GeoMesh.h"
#include "tpanic.h"
#include "GeoElementSide.h"
#include "Poisson.h"

GeoMesh::GeoMesh(){
   
}

GeoMesh::GeoMesh(const GeoMesh & cp){
    this->operator =(cp);
}

GeoMesh &GeoMesh::operator=(const GeoMesh &cp){
    
    int nel = cp.Elements.size();
    int nnodes = cp.Nodes.size();
    
    this->Nodes.resize(nnodes);
    for(int inodes = 0; inodes < nnodes; inodes++)
    {
        this->Nodes[inodes] = cp.Nodes[inodes];
    }
    
    this->Elements.resize(nel);
    for(int iel = 0; iel < nel ; iel++)
    {
        if (cp.Elements[iel])
        {
            this->Elements[iel] = cp.Elements[iel];
        }
        else
        {
            this->Elements[iel] = NULL;
        }
    }
    return *this;
}

void GeoMesh::SetNumNodes(int nnodes){
    Nodes.resize(nnodes);
}

void GeoMesh::SetNumElements(int numelements){
    Elements.resize(numelements);
}

int GeoMesh::NumNodes(){
    return Nodes.size();
}

int GeoMesh::NumElements(){
    return Elements.size();
}

GeoNode &GeoMesh::Node(int node){
    return Nodes[node];
}

void GeoMesh::SetElement(int elindex, GeoElement *gel){
    Elements[elindex]=gel;
}

GeoElement *GeoMesh::Element(int elindex){
    return Elements[elindex];
}

void GeoMesh::BuildConnectivity(){
    
    
            int NumNds = NumNodes();
            int NumEl =NumElements();
            for(int i= 0; i<NumElements(); i++){
                GeoElement *geoel = Element(i);
                geoel->SetIndex(i);
            }
    
    
            Matrix Conects(0,3);
            std::vector<Matrix> VecConec(NumNds);
            //verifica que elementos comparten el mismo nodo y los almacena en Neig
            for(int i=0; i< NumNds; i++){
                std::vector<int> Neig;
                for(int j =0; j< NumEl;j++){
                    GeoElement *element = Element(j);
                    int NodsEl= element->NNodes();
                    for(int k=0; k<= NodsEl; k++){
                        int64_t nodindex = element->NodeIndex(k);
                        if(nodindex==i){
                            //        std::cout<<"Elemento: "<<j<<"Index :"<<nodindex<<"\n";
                            Neig.push_back(j);
    
                        }
                    }
                }
    
                //guarda la informacion de Neig en una matriz
                int numNeig = Neig.size();
                Conects.Resize(numNeig, 3);
                for(int w=0; w<numNeig; w++){
                    GeoElement *element = Element(Neig[w]);
                    int elSides= element->NSides();
                    for(int l=0; l<elSides; l++){
                        int64_t nodindex = element->NodeIndex(l);
                        if(i==nodindex){
//                            std::cout<<"Element: "<< Neig[w] <<" side: "<<l<<" node:"<< nodindex <<"\n";
                            Conects(w,0)=Neig[w];
                            Conects(w,1)=l;
                            Conects(w,2)=nodindex;
    
                        }
                    }
                }
//                std::cout<<"Elementos asociados al nodo: "<<i<<"\n";
//                Conects.Print();
                VecConec[i]=Conects;
            }
//        VecConec[0].Print();
    
            // Conectividad por ramales por nodos
            Matrix NodoConect;
            int NumconectsNo;
            for(int no=0; no<NumNds; no++){
                NodoConect = VecConec[no];
                NumconectsNo = NodoConect.Rows();
                for(int i=0; i<(NumconectsNo-1); i++ ){
                    int idEl = NodoConect(i+1,0);
                    int idElAn = NodoConect(i,0);
                    int side = NodoConect(i+1,1);
                    int sideAn = NodoConect(i,1);
                    GeoElement *geoel = Element(idEl);
                    GeoElement *geoelAn = Element(idElAn);
                    GeoElementSide geoside(geoel,side);
                    geoelAn->SetNeighbour(sideAn, geoside);
                }
                int idEl = NodoConect(0,0);
                int idElAn = NodoConect(NumconectsNo -1, 0) ;
                int side = NodoConect(0,1);
                int sideAn = NodoConect(NumconectsNo -1, 1);
                GeoElement *geoel = Element(idEl);
                GeoElement *geoelAn = Element(idElAn);
                GeoElementSide geoside(geoel,side);
                geoelAn->SetNeighbour(sideAn, geoside);
            }
    
    
    for(int no=0; no<NumNds; no++){
         NodoConect = VecConec[no];
         int numNeih = NodoConect.Rows();
            for(int i=0; i< numNeih; i++){
                GeoElement *geoelAn = Element( NodoConect(i,0));
                VecInt VecNode1;
                geoelAn->GetNodes(VecNode1);
                for(int j=0; j< numNeih; j++){
                    if(i!=j){
                    GeoElement *geoel = Element( NodoConect(j,0));
                    VecInt VecNode2;
                    geoel->GetNodes(VecNode2);
                
                    VecInt NodIntersection;
                    std::sort(VecNode1.begin(),VecNode1.end());
                    std::sort(VecNode2.begin(),VecNode2.end());
                    std::set_intersection(VecNode1.begin(), VecNode1.end(), VecNode2.begin(), VecNode2.end(),std::back_inserter(NodIntersection));
                        if (NodIntersection.size()==2){
                            int sideAn = geoelAn->WhichSide(NodIntersection);
                            int side =   geoel->WhichSide(NodIntersection);
                            GeoElementSide GeoSideAn(geoelAn,sideAn);
                            GeoElementSide GeoSide(geoel,side);
                            geoelAn->SetNeighbour(sideAn, GeoSide);
                           // geoel->SetNeighbour(side, GeoSideAn);
                        
                        
                        }
                    }
                }
        }
    }
    
    
    
    
    
    
    
    
   
    
   
    
 }

void GeoMesh::Print(std::ostream &out){
    int i;
    out << "PRINTING GEOMESH..." << std::endl;
    out << "Numero de Nodos = " << Nodes.size() << std::endl;
    for (i=0;i<Nodes.size();i++)
    {
        out << "Index " << i << ' ';
        Nodes[i].Print(out);
    }
    out <<"\n" <<"Numero de Elementos = " << Elements.size() << std::endl;
    for (i=0; i<this->Elements.size(); i++)
    {
        out << "Element Index " << i << ' ';
        Elements[i]->Print(out);
    }
    
}
