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
    fSide = 0;
    fElement = 0;
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


// Fill in the data structure for the neighbouring information
void SetNeighbour(const GeoElementSide &neighbour){
    
}

// Verifiy if an element is a neighbour
bool IsNeighbour(const GeoElementSide &candidate){
    
}

// Define elements neighbourhood
void IsertConnectivity(GeoElementSide &connectivity){
    
}

// Vector with all Neighbours
void AllNeighbours(std::vector<GeoElementSide> &allneigh){
    
}

// Compute all corner neighbours
void ComputeNeighbours(std::vector<GeoElementSide> &neighbour){
    
}
