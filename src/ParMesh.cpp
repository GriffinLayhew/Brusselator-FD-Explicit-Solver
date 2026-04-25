#include "ParMesh.hpp"

// Constructor Definition
ParMesh::ParMesh(int rank, int size, int N, double Domain, int numColumns)
{
    myRank      = rank;
    globalSize  = size;

    // int test = ensureGoodDecomp(globalSize, numColumns);     NEED TO FIX

    // Determine the Amount of Data my Chunk is Repsonsible for
    int numRows = size / numColumns;
    int myNx = N / numColumns;
    int myNy = N / numRows;

    int totalGridPoints = myNx * myNy;

    rowColIndex myIndicies = getMyIndexes(myRank, globalSize, numColumns);

    // Left Neighboor
    rowColIndex leftIndicies;
    leftIndicies.Row = myIndicies.Row;
    leftIndicies.Column = myIndicies.Column - 1;

    int leftRank = getRankfromIndicies(leftIndicies, numColumns, numRows);
    neighboorRanks.push_back(leftRank);

    // Right Neighboor
    rowColIndex rightIndicies;
    rightIndicies.Row = myIndicies.Row;
    rightIndicies.Column = myIndicies.Column + 1;

    int rightRank = getRankfromIndicies(rightIndicies, numColumns, numRows);
    neighboorRanks.push_back(rightRank);

    // Top Neighboor
    rowColIndex topIndicies;
    topIndicies.Row = myIndicies.Row - 1;
    topIndicies.Column = myIndicies.Column;

    int topRank = getRankfromIndicies(topIndicies, numColumns, numRows);
    neighboorRanks.push_back(topRank);

    // Bottom Neighboor
    rowColIndex bottomIndicies;
    bottomIndicies.Row = myIndicies.Row + 1;
    bottomIndicies.Column = myIndicies.Column;

    int bottomRank = getRankfromIndicies(bottomIndicies, numColumns, numRows);
    neighboorRanks.push_back(bottomRank);

    // Fill my Chunk with Zeros with Halos
    int totalWithHalos = (myNy + 2) * (myNx + 2);
    for(int i = 0; i < totalWithHalos; i++)
        data.push_back(rank);
    

}

// Non-Modifying Definitions
void ParMesh::display(int numColumns, int N) const
{
    // Determine the amount of data my chunk is responsible for
    int numRows = this->getSize() / numColumns;
    int myNx = N / numColumns;
    int myNy = N / numRows;

    int stride = myNx + 2;   // row width including left/right halos

    // Formatting controls
    const int precision = 3;
    const int width = 8;

    // ---------------- DISPLAY MESH INFORMATION ----------------
    std::cout << "\n\nRank = " << myRank
              << "   --------------   "
              << myNx << " x " << myNy << "\n\n";

    std::cout << std::fixed << std::setprecision(precision);

    // ---------------- TOP HALO (no corners) ----------------
    std::cout << std::setw(width) << " " << std::setw(width) << " ";
    for (int j = 1; j <= myNx; j++)
        std::cout << std::setw(width) << data[j];
    std::cout << '\n';

    // ---------------- TOP BARRIER ----------------
    std::cout << std::setw(width) << " ";
    for (int j = 0; j < myNx + 2; j++)
        std::cout << std::setw(width) << "-----";
    std::cout << '\n';

    // ---------------- INTERIOR ROWS WITH SIDE HALOS ----------------
    for (int i = 1; i <= myNy; i++)
    {
        int rowStart = i * stride;

        // left halo
        std::cout << std::setw(width) << data[rowStart];

        // left separator
        std::cout << std::setw(width) << "|";

        // interior
        for (int j = 1; j <= myNx; j++)
            std::cout << std::setw(width) << data[rowStart + j];

        // right separator and halo
        std::cout << std::setw(width) << "|"
                  << std::setw(width) << data[rowStart + myNx + 1]
                  << '\n';
    }

    // ---------------- BOTTOM BARRIER ----------------
    std::cout << std::setw(width) << " ";
    for (int j = 0; j < myNx + 2; j++)
        std::cout << std::setw(width) << "-----";
    std::cout << '\n';

    // ---------------- BOTTOM HALO (no corners) ----------------
    int bottomStart = (myNy + 1) * stride;
    std::cout << std::setw(width) << " " << std::setw(width) << " ";
    for (int j = 1; j <= myNx; j++)
        std::cout << std::setw(width) << data[bottomStart + j];
    std::cout << '\n';
}

double ParMesh::getData(int i, int j, int myNx, int myNy) const
{
    int stride = myNx + 2;
    int index = i * stride + j;
    return data[index];
}


// Modifying Member Functions
void ParMesh::haloExchange(int numColumns, int N)
{

    // Get my Amount to Send and Recieve
    int numRows = this->getSize() / numColumns;
    int myNx = N / numColumns;
    int myNy = N / numRows;

    int stride = myNx + 2;     // row width including left/right halos

    // Retrieve My Neighboors
    int rankLeft    = neighboorRanks[0];
    int rankRight   = neighboorRanks[1];
    int rankTop     = neighboorRanks[2];
    int rankBottom  = neighboorRanks[3];

    // Create the Send and Recieve Buffers for Each
    double* dataLeft_Send = new double[myNy];
    double* dataLeft_Recv = new double[myNy];

    double* dataRight_Send = new double[myNy];
    double* dataRight_Recv = new double[myNy];

    double* dataTop_Send = new double[myNx];
    double* dataTop_Recv = new double[myNx];


    double* dataBottom_Send = new double[myNx];
    double* dataBottom_Recv = new double[myNx];

    // Fill up Send Buffers
        // ------TOP
        for (int j = 1; j <= myNx; j++)
            dataTop_Send[j-1] = data[j];
        
        // ------BOTTOM
        int bottomStart = (myNy) * stride;
        int count = 0;
        for (int j = 1; j <= myNx; j++)
        {
            dataBottom_Send[count] = data[bottomStart + j];
            count++;
        }
            

        //  -----LEFT AND RIGHT
        for (int i = 1; i <= myNy; i++)
        {
            int rowStart = i * stride;

            // left halo
            dataLeft_Send[i - 1] = data[rowStart + 1];

            // right halo
            dataRight_Send[i - 1] = data[rowStart + myNx];

        }


        //-----------------HORIZONTAL SEND AND RECIEVE-------------------------//
    MPI_Request yReq[4];
    MPI_Status  ymyStatus[4];

    for (int i = 0; i < 4; i++)
        yReq[i] = MPI_REQUEST_NULL;

    // Exchange Left
    if (rankLeft != -1)
    {
        MPI_Irecv(dataLeft_Recv, myNy, MPI_DOUBLE, rankLeft, 100, MPI_COMM_WORLD, &yReq[0]);
        MPI_Isend(dataLeft_Send, myNy, MPI_DOUBLE, rankLeft, 101, MPI_COMM_WORLD, &yReq[1]);
    }
    else
    {
        int count = 0;
        for (int i = 1; i <= myNy; i++)
        {
            int rowStart = i * stride;
            dataLeft_Recv[count] = data[rowStart + 1];
            count++;
        }
    }

    // Exchange Right
    if (rankRight != -1)
    {
        MPI_Irecv(dataRight_Recv, myNy, MPI_DOUBLE, rankRight, 101, MPI_COMM_WORLD, &yReq[2]);
        MPI_Isend(dataRight_Send, myNy, MPI_DOUBLE, rankRight, 100, MPI_COMM_WORLD, &yReq[3]);
    }
    else
    {
        int count = 0;
        for (int i = 1; i <= myNy; i++)
        {
            int rowStart = i * stride;
            dataRight_Recv[count] = data[rowStart + myNx];
            count++;
        }
    }

    MPI_Waitall(4, yReq, ymyStatus);

        //-----------------VERTICAL SEND AND RECIEVE-------------------------//
    MPI_Request xReq[4];
    MPI_Status  xmyStatus[4];

    for (int i = 0; i < 4; i++)
        xReq[i] = MPI_REQUEST_NULL;

    // Exchange Top
    if (rankTop != -1)
    {
        MPI_Irecv(dataTop_Recv, myNx, MPI_DOUBLE, rankTop, 200, MPI_COMM_WORLD, &xReq[0]);
        MPI_Isend(dataTop_Send, myNx, MPI_DOUBLE, rankTop, 201, MPI_COMM_WORLD, &xReq[1]);
    }
    else
    {
        for (int j = 1; j <= myNx; j++)
            dataTop_Recv[j-1] = data[j];
    }

    // Exchange Bottom
    if (rankBottom != -1)
    {
        MPI_Irecv(dataBottom_Recv, myNx, MPI_DOUBLE, rankBottom, 201, MPI_COMM_WORLD, &xReq[2]);
        MPI_Isend(dataBottom_Send, myNx, MPI_DOUBLE, rankBottom, 200, MPI_COMM_WORLD, &xReq[3]);
    }
    else
    {
        int bottomStart = (myNy) * stride;
        int count = 0;
        for (int j = 1; j <= myNx; j++)
        {
            dataBottom_Recv[count] = data[bottomStart + j];
            count++;
        }
    }

    MPI_Waitall(4, xReq, xmyStatus);

    // Fill Exhanged Data into ParMesh's 'data'

        // Fill Left & Right

    for (int i = 1; i <= myNy; i++)
    {
        int rowStart = i * stride;

        // left halo
        data[rowStart] = dataLeft_Recv[i - 1];

        // right halo
        data[rowStart + myNx + 1] = dataRight_Recv[i - 1];
    }


        
    for (int j = 1; j <= myNx; j++)
    {
        int bottomHaloStart = (myNy + 1) * stride;

        // top halo
        data[j] = dataTop_Recv[j - 1];

        // bottom halo
        data[bottomHaloStart + j] = dataBottom_Recv[j - 1];
    }
    
    


    // Delete and Clear up Memory
    delete[] dataLeft_Send;
    delete[] dataLeft_Recv;

    delete[] dataRight_Send;
    delete[] dataRight_Recv;

    delete[] dataTop_Send;
    delete[] dataTop_Recv;

    delete[] dataBottom_Send;
    delete[] dataBottom_Recv;
}

void ParMesh::updateData(int i, int j, int myNx, int myNy, double newData)
{
    int stride = myNx + 2;
    int index = i * stride + j;
    data[index] = newData;
}


// Helper Definitions
rowColIndex getMyIndexes(int myRank, int globalSize, int numColumns)
{
    rowColIndex myIndexes;

    myIndexes.Column    = myRank % numColumns;
    myIndexes.Row       = myRank / numColumns;

    return myIndexes;
}

int getRankfromIndicies(rowColIndex Indicies, int numColumns, int numRows)
{
    if (Indicies.Row >= numRows || Indicies.Row < 0 || Indicies.Column >= numColumns || Indicies.Column < 0)
        return -1;
    
    return Indicies.Row*numColumns + Indicies.Column;
}

int ensureGoodDecomp(int &globalSize, int &numColumns)
{
    // Basic sanity
    if (globalSize <= 0 || numColumns <= 0)
        return -1;

    // numColumns cannot exceed number of ranks
    if (numColumns > globalSize)
        numColumns = globalSize;

    // If already valid, nothing to do
    if (globalSize % numColumns == 0)
        return 0;

    // Find nearest divisor of globalSize to numColumns
    int bestCols = 1;
    int bestDist = std::abs(numColumns - 1);

    for (int c = 1; c <= globalSize; ++c)
    {
        if (globalSize % c == 0)
        {
            int dist = std::abs(numColumns - c);

            // choose closest divisor; if tied, prefer smaller
            if (dist < bestDist || (dist == bestDist && c < bestCols))
            {
                bestCols = c;
                bestDist = dist;
            }
        }
    }

    numColumns = bestCols;
    return 0;
}

