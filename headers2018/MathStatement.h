//
//  MathStatement.h
//  FemCourse
//
//  Created by Philippe Devloo on 24/04/18.
//

#ifndef MathStatement_h
#define MathStatement_h

#include "DataTypes.h"
#include "IntPointData.h"

class MathStatement
{
    // Math statement ID
    int matid = 0;
    
    // Number of state variable
    int nstate = 0;
    
public:
    
    // Constructor of MathStatement
    MathStatement();
    
    // Copy constructor of MathStatement
    MathStatement(const MathStatement &copy);
    
    // Operator of copy
    MathStatement &operator=(const MathStatement &copy);
    
    // Destructor of MathStatement
    virtual ~MathStatement();
    
    // Method for creating a copy of the element
    virtual MathStatement *Clone() const = 0;
    
    // Return the number of state variables
    virtual int NState() const = 0;
    
    // Method to implement integral over element's volume
    virtual void Contribute(IntPointData &integrationpointdata, Matrix &EK, Matrix &EF) const = 0;
    
    // Method to compute the contribution to the error norm
    virtual void ContributeError(IntPointData &integrationpointdata, std::function<void(const VecDouble &co, VecDouble &sol, Matrix &dsol)> &exact)=0;
    
};
#endif /* MathStatement_h */
