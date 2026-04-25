#pragma once
#include "Utilities.hpp"

typedef struct 
{
    int Row;
    int Column;

} rowColIndex;



class ParMesh
{
    public:
    // Constructor
    ParMesh(int rank, int size, int N, double Domain, int numColumns);

    // Non-Modifying Member Functions
    std::vector<int>     getNeighboorRanks() const {return neighboorRanks;}
    int                  getRank()           const {return myRank;}
    int                  getSize()           const {return globalSize;}
    void                 display(int, int)   const; 
    double               getData(int j, int k, int, int)      const;
    int                  getNumEntries()     const  {return data.size();}
    std::vector<double>& getRawData()                 {return data;}

    // Modifying Member Functions
    void                haloExchange(int, int);
    void                updateData(int, int, int, int, double);

    private:
        int                     myRank;
        int                     globalSize;
        std::vector<int>        neighboorRanks; // [L, R, T, B]
        std::vector<double>     data;

};

// Helper Functions

rowColIndex getMyIndexes(int myRank, int globalSize, int numColumns);
int         getRankfromIndicies(rowColIndex Indicies, int numColumns, int numRows);
int         ensureGoodDecomp(int &globalSize, int &numColumns);
