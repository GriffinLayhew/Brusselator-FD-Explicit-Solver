#include "State.hpp"

// Constructor Definition
State::State(int rank, int size, int N, double Domain, int numColumns, std::vector<std::string> varNames)
{
    for (int i=0; i < varNames.size(); i++)
    {
        Field newField(rank, size, N, Domain, numColumns, varNames[i]);
        scalarVariables.push_back(newField);
    }
    
}

// Non-Modifying Definitions


// Modifying Member Functions



// Helper Definitions

