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
                    std::cout<<"Element: "<< Neig[w] <<" side: "<<l<<" node:"<< nodindex <<"\n";
                    Conects(w,0)=Neig[w];
                    Conects(w,1)=l;
                    Conects(w,2)=nodindex;
                 
                }
            }
        }
        VecConec[i]=Conects;
    }
    VecConec[0].Print();
    
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
    //Conectividad 1D
    for(int i=0; i< NumEl; i++){
        GeoElement *geoel = Element(i);
        int nNodes = geoel->NNodes();
        
        
        
        
    }
    
    // Conectividad 2D.
    for(int i=0; i<= NumEl; i++){
        
    }
    
    
    
    

////Vector de Vecinos de todos los elementos
//    std::vector< std::vector < std::vector <GeoElementSide> >  > vecinos(NumEl);
//
//    for(int i=0; i< NumEl; i++){
//        GeoElement *element = Element(i);
//        int nsidesEli = element->NSides();
//    vecinos[i].resize(nsidesEli);
//    }
//
//  //setar vecinos por informacion de los nodos
//    for(int i=0; i< NumNds; i++){
//        Matrix vecNod = VecConec[i];
//        int tam = vecNod.Rows();
//        for (int j=0; j<tam;j++) {
//            int lado = vecNod(j,1);
//            for (int k=0; k<tam; k++){
//                if(j!=k){
//                    int ubElement = vecNod(k,0);
//                    GeoElement *element = Element(ubElement);
//                    int ladoVeci = vecNod(k,1);
//                    GeoElementSide neigh(element,ladoVeci);
//                    int elementAn = vecNod(j,0);
////                    std::cout<<"vecinos del elemento: "<<elementAn<< " por el lado "<<lado<<"\n";
////                    std::cout<<" "<<ubElement<< "/"<<ladoVeci<<"\n";
//                    vecinos[elementAn][lado].push_back(neigh);
//                }
//            }
//        }
//
//    }
////conectividad por lados (elementos 2D con el elemento lineal)
//    for(int i=0; i<NumEl; i++){
//        GeoElement *elementi = Element(i);
//        int Eltypei = elementi->Type();
//        if(Eltypei==1){
//            int no1=elementi->NodeIndex(0);
//            int no2=elementi->NodeIndex(1);
//
//            for(int j=0; j<NumEl;j++){
//                GeoElement *elementj = Element(j);
//                int Eltypej = elementj->Type();
//                if(Eltypej!=1){
//                    int cornersj= elementj->NCornerNodes();
//
//                    //HERE1
////Cierra el elemento 2D
//                    int noj1=elementj->NodeIndex(cornersj-1);
//                    int noj2=elementj->NodeIndex(0);
////side es el ultimo corner del elemento 2D lo une con la primer corner de el mismo
//                    int side= elementj->NSides()-2;
//                    GeoElementSide neigh(elementj,side);
//                    GeoElementSide neigh2(elementi,2);
//                    if(no1==noj1 ||  no1==noj2){
//                        if(no2==noj2 ||  no2==noj1 ){
//                            vecinos[i][2].push_back(neigh);
//                            vecinos[j][side].push_back(neigh2);
////                          std::cout<<"vecino del elemento: "<<i <<" por el lado: "<<2<<"\n";
////                            std::cout<< j<<"/"<<side<<"\n";
////
////                            std::cout<<"vecino del elemento: "<<j <<" por el lado: "<<side<<"\n";
////                            std::cout<< i<<"/"<<"2"<<"\n";
//                        }
//                    }
//
//
//                    //HERE2
//
////CALCULA LADO POR LADO MENOS EL LADO PARA CERRAR
//                    for(int k=0; k<cornersj-1;k++){
//                        int noj1=elementj->NodeIndex(k);
//                        int noj2=elementj->NodeIndex(k+1);
//                        int side= k + 4;
//                        GeoElementSide neigh(elementj,side);
//                        GeoElementSide neigh2(elementi,2);
//
//                        if(no1==noj1 ||  no1==noj2){
//                            if(no2==noj2 ||  no2==noj1 ){
//
//                                vecinos[i][2].push_back(neigh);
//                                vecinos[j][side].push_back(neigh2);
//                                std::cout<<"vecino del elemento: "<<i <<" por el lado: "<<2<<"\n";
//                                std::cout<< j<<"/"<<side<<"\n";
//
//                                std::cout<<"vecino del elemento: "<<j <<" por el lado: "<<side<<"\n";
//                                std::cout<< i<<"/"<<"2"<<"\n";
//                            }
//                            int numVecinos =vecinos[i][2].size();
//                            if(numVecinos==2){
//                                for(int z=0; z < NumEl; z++){
//                                    GeoElement *elementz = Element(z);
////                                    GeoElement *elementv1 = vecinos[i][2][0].fElement;
//                                    int side1 =vecinos[i][2][0].fSide;
////                                    GeoElement *elementv2 = vecinos[i][2][1].fElement;
//                                     int side2 =vecinos[i][2][1].fSide;
//
//
//                                    if(elementz == vecinos[i][2][0].fElement){
//                                    vecinos[z][side1].push_back(vecinos[i][2][1]);
//                                    }
//                                    if(elementz == vecinos[i][2][1].fElement){
//                                        vecinos[z][side2].push_back(vecinos[i][2][0]);
//                                    }
//
//                                    }
//                            }
//                        }
//
//                    }
//
//                }
//            }
//
//        }
//       // std::cout;
//
//    }
//
//
//
    
    //hasta aqui comentamos
    
//
//    //solo imprime las conectividades
//    for(int i=0;i<NumEl;i++){
//        std::cout<<"Vecinos del elemento: "<<i<<"\n";
//        GeoElement *elementi = Element(i);
//        int numsides = elementi->NSides();
//        for (int j=0; j<numsides; j++) {
//            std::cout<<"side "<<j<<" : ";
//            int nVecinos = vecinos[i][j].size();
//            if(nVecinos==0){
//                std::cout<<" No neighbour";
//            }
//            for(int k=0; k< nVecinos; k++){
//               GeoElement *elementk = vecinos[i][j][k].fElement;
//                int side = vecinos[i][j][k].fSide;
//                for(int w=0; w< NumEl;w++){
//                    if(Element(w)==elementk){
//                        std::cout<<" "<<w<<"/"<<side;
//                    }
//                }
//            }
//            std::cout<<"\n";
//        }
//
//    }
    
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
