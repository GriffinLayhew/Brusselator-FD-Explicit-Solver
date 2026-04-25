#include "Output.hpp"

// Constructor Definition



// Non-Modifying Definitions
void Output::printToFile()
{
    int rank = myState.getField(1).getCurrent().getRank();
    int size = myState.getField(1).getCurrent().getSize();

    int NumColumns = myState.getField(0).getCols();
    int N          = myState.getField(0).getN();

    int numRows = size / NumColumns;
    int myNx = N / NumColumns;
    int myNy = N / numRows;

    int numEntriesLocal = myNx * myNy;
    double* allEntriesLocal = new double[numEntriesLocal];

    int idx = 0;
    for (int i = 1; i <= myNy; i++)
    {
        for (int j = 1; j <= myNx; j++)
        {
            allEntriesLocal[idx++] =
                myState.getField(0).getCurrent().getData(i, j, myNx, myNy);
        }
    }

    std::vector<int> all_sizes(size, numEntriesLocal);
    std::vector<int> displacements(size);

    for (int r = 0; r < size; r++)
        displacements[r] = r * numEntriesLocal;

    double* gathered = nullptr;

    if (rank == 0)
        gathered = new double[N * N];

    MPI_Gatherv(allEntriesLocal,
                numEntriesLocal,
                MPI_DOUBLE,
                gathered,
                all_sizes.data(),
                displacements.data(),
                MPI_DOUBLE,
                0,
                MPI_COMM_WORLD);

    if (rank == 0)
    {
        std::vector<double> globalGrid(N * N, 0.0);

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
                    int localIndex  = base + i * myNx + j;

                    globalGrid[globalIndex] = gathered[localIndex];
                }
            }
        }

        std::ofstream MyFile(fileName, std::ios::app);

        if (MyFile.is_open())
        {
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    MyFile << std::setw(10)
                           << std::setprecision(5)
                           << std::fixed
                           << globalGrid[i * N + j];
                }
                MyFile << "\n";
            }

            MyFile << "\n\n\n";
            MyFile.close();
        }
        else
        {
            std::cerr << "Unable to open file\n";
        }

        delete[] gathered;
    }

    delete[] allEntriesLocal;
}

// Modifying Member Functions



// Helper Definitions

