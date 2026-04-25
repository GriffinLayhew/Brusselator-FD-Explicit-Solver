#pragma once

#include "Utilities.hpp"
#include "Field.hpp"

class State
{
    public:
    // Constructor
        State(int rank, int size, int N, double Domain, int numColumns, std::vector<std::string> varNames);


    // Non-Modifying Member Functions
        Field& getField(int i)           {return scalarVariables[i];}
        int    getNumdFields()   const   {return scalarVariables.size();}

    // Modifying Member Functions

    private:
        std::vector<Field> scalarVariables;

};

// Helper Functions


