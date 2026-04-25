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
State State::operator+(State other)
{
    State result = *this; // copy structure + data

    int numFields = this->getNumdFields();

    for (int k = 0; k < numFields; k++)
    {
        const std::vector<double>& dataA =
            this->getField(k).getCurrent().getRawData();

        const std::vector<double>& dataB =
            other.getField(k).getCurrent().getRawData();

        std::vector<double>& dataR =
            result.getField(k).getCurrent().getRawData();

        int size = dataA.size();

        for (int i = 0; i < size; i++)
        {
            dataR[i] = dataA[i] + dataB[i];
        }
    }

    return result;
}


// Helper Definitions

