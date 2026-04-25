#include "Field.hpp"

// Constructor Definition
Field::Field(int rank, int size, int N, double Domain, int numColumns, std::string myFieldName) : currentIteration(rank, size, N, Domain, numColumns), 
                                                                                                  nextIteration(rank, size, N, Domain, numColumns),
                                                                                                  Name(myFieldName)
{
    NumCols = numColumns;
    Ntotal = N;
    DomainSize = Domain;
}

// Non-Modifying Definitions


// Modifying Member Functions



// Helper Definitions

