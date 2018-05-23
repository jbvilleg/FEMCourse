//
//  MathStatement.h
//  FemCourse
//
//  Created by Philippe Devloo on 24/04/18.
//

#include "MathStatement.h"


MathStatement::MathStatement(){
    
}

MathStatement::MathStatement(const MathStatement &copy){
    nstate=copy.nstate;
    matid=copy.matid;
}

MathStatement &MathStatement::operator=(const MathStatement &copy){
    nstate=copy.nstate;
    matid=copy.matid;
    return *this;
}

MathStatement::~MathStatement(){
    
}

