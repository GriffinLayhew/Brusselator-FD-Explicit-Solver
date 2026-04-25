#include "Output.hpp"

// Constructor Definition



// Non-Modifying Definitions
void Output::printToFile()
{
    std::ofstream MyFile(fileName, std::ios::app);

    int numFields = myState.getNumdFields();
    int rank = myState.getField(1).getCurrent().getRank();
    int size = myState.getField(1).getCurrent().getSize();

    ParMesh currentU = myState.getField(0).getCurrent();
    int NumColumns   = myState.getField(0).getCols();
    int N            = myState.getField(0).getN();

    int numRows = size / NumColumns;
    int myNx = N / NumColumns;
    int myNy = N / numRows;
    int stride = myNx + 2;

    int numEntriesLocal = myNx * myNy;
    double* allEntriesLocal = new double[numEntriesLocal];

    // ---------------- PACK LOCAL DATA (NO HALOS) ----------------
    int idx = 0;
    for (int i = 1; i <= myNy; i++)
    {
        for (int j = 1; j <= myNx; j++)
        {
            allEntriesLocal[idx++] =
                myState.getField(1).getCurrent().getData(i, j, myNx, myNy);
        }
    }

    // ---------------- BUILD GATHERV ARRAYS ----------------
    std::vector<int> all_sizes(size, numEntriesLocal);
    std::vector<int> displacements(size);

    for (int r = 0; r < size; r++)
        displacements[r] = r * numEntriesLocal;

    int numEntriesGlobal = N * N;
    double* allEntriesPrintArray = nullptr;

    if (rank == 0)
        allEntriesPrintArray = new double[numEntriesGlobal];

    // ---------------- GATHER ----------------
    MPI_Gatherv(allEntriesLocal,
                numEntriesLocal,
                MPI_DOUBLE,
                allEntriesPrintArray,
                all_sizes.data(),
                displacements.data(),
                MPI_DOUBLE,
                0,
                MPI_COMM_WORLD);

    // ---------------- WRITE FILE ----------------
    if (rank == 0)
    {
        if (MyFile.is_open())
        {
            // Reconstruct global 2D layout
            for (int r = 0; r < size; r++)
            {
                rowColIndex idxRC = getMyIndexes(r, size, NumColumns);

                int rowOffset = idxRC.Row * myNy;
                int colOffset = idxRC.Column * myNx;

                int base = r * numEntriesLocal;

                for (int i = 0; i < myNy; i++)
                {
                    for (int j = 0; j < myNx; j++)
                    {
                        int globalRow = rowOffset + i;
                        int globalCol = colOffset + j;

                        int globalIndex = globalRow * N + globalCol;

                        allEntriesPrintArray[globalIndex] =
                            allEntriesPrintArray[base + i * myNx + j];
                    }
                }
            }

            // Print as 2D grid
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    MyFile << std::setw(10)
                           << std::setprecision(5)
                           << std::fixed
                           << allEntriesPrintArray[i * N + j];
                }
                MyFile << "\n";
            }

            
        }
        else
        {
            std::cerr << "Unable to open file\n";
        }

        std::cout << "Successfully written to file.\n";
        MyFile << std::endl << std::endl << std::endl;
        MyFile.close();
        
    }

    // ---------------- CLEANUP ----------------
    delete[] allEntriesLocal;
    if (rank == 0)
        delete[] allEntriesPrintArray;
}

// Modifying Member Functions



// Helper Definitions

