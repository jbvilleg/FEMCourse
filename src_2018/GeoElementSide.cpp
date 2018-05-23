//
//  Geom1d.h
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//


#include "Topology1d.h"
#include "Geom1d.h"
#include "GeoElementSide.h"


GeoElementSide::GeoElementSide(){
//    fSide = -1;
}



GeoElementSide::GeoElementSide(const GeoElementSide &copy){
    fElement=copy.Element();
    fSide=copy.fSide;
}

GeoElementSide &GeoElementSide::operator=(const GeoElementSide &copy){
    this->fElement=copy.fElement;
    fSide = copy.fSide;
    return *this;
}



GeoElementSide GeoElementSide::Neighbour() const{
    return fElement ? fElement->Neighbour(fSide) : GeoElementSide();
}


/** @brief Fill in the data structure for the neighbouring information*/
void GeoElementSide::SetNeighbour( GeoElementSide &neighbour){
    fElement->SetNeighbour(fSide,neighbour);
     
}
