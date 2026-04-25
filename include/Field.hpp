#pragma once

#include "Utilities.hpp"
#include "ParMesh.hpp"

class Field
{
    public:
    // Constructor
    Field(int rank, int size, int N, double Domain, int numColumns, std::string myFieldName);

    // Non-Modifying Member Functions
    ParMesh&       getCurrent()       { return currentIteration; }
    ParMesh&       getNext()          { return nextIteration; }

    int     getCols()       {return NumCols;}
    int     getN()          {return Ntotal;}
    double  getDomain()     {return DomainSize;}


    // Modifying Member Functions

    private:
        ParMesh         currentIteration;
        ParMesh         nextIteration;
        std::string     Name;
        int             NumCols;
        int             Ntotal;
        double          DomainSize;

};

// Helper Functions

