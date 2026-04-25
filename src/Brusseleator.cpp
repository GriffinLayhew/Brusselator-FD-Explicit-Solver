#include "Brusselator.hpp"

// Constructor Definition
Brusselator::Brusselator(Parameters params, State& myState) : myParams{params}, myVars{myState}
{

}

// Non-Modifying Definitions
void Brusselator::display()
{
    ParMesh currentU = myVars.getField(0).getCurrent();
    ParMesh currentV = myVars.getField(1).getCurrent();

    std::cout << "===================== U =====================" << std::endl;
    currentU.display(myVars.getField(0).getCols(), myVars.getField(0).getN());
    std::cout << "===================== V =====================" << std::endl;
    currentV.display(myVars.getField(0).getCols(), myVars.getField(0).getN());

}


// Modifying Member Functions
void Brusselator::forwardEuler(double coeff)
{
    ParMesh& currentU = myVars.getField(0).getCurrent();
    ParMesh& nextU    = myVars.getField(0).getNext();

    ParMesh& currentV = myVars.getField(1).getCurrent();
    ParMesh& nextV    = myVars.getField(1).getNext();

    int NumColumns = myVars.getField(0).getCols();
    int N          = myVars.getField(0).getN();

    int numRows = currentU.getSize() / NumColumns;
    int myNx = N / NumColumns;
    int myNy = N / numRows;
    int stride = myNx + 2;     // row width including left/right halos

    double dx = myVars.getField(0).getDomain() / N;
    double dy = myVars.getField(0).getDomain() / N;



    // ---------------- INTERIOR ROWS WITH SIDE HALOS ----------------
    for (int i = 1; i <= myNy; i++)
    {
        for (int j = 1; j <= myNx; j++)
        {
            double Diff_U = myParams.Du;
            double Diff_V = myParams.Dv;
            double A = myParams.A;
            double B = myParams.B;

            double UnewDiff = Diff_U * ( (currentU.getData(i+1,j,myNx,myNy) - 2 * currentU.getData(i,j,myNx,myNy) + currentU.getData(i-1,j,myNx,myNy)) / (dx*dx) + (currentU.getData(i,j+1,myNx,myNy) - 2 * currentU.getData(i,j,myNx,myNy) + currentU.getData(i,j-1,myNx,myNy)) / (dy*dy));
            double VnewDiff = Diff_V * ( (currentV.getData(i+1,j,myNx,myNy) - 2 * currentV.getData(i,j,myNx,myNy) + currentV.getData(i-1,j,myNx,myNy)) / (dx*dx) + (currentV.getData(i,j+1,myNx,myNy) - 2 * currentV.getData(i,j,myNx,myNy) + currentV.getData(i,j-1,myNx,myNy)) / (dy*dy));

            double updateU = currentU.getData(i,j,myNx,myNy) +  coeff * (UnewDiff + A + currentU.getData(i,j,myNx,myNy) * currentU.getData(i,j,myNx,myNy) * currentV.getData(i,j,myNx,myNy) - (B+1)*currentU.getData(i,j,myNx,myNy));
            nextU.updateData(i, j, myNx, myNy, updateU); 
        
            double updateV = currentV.getData(i,j,myNx,myNy) + coeff * (VnewDiff + B*currentU.getData(i,j,myNx,myNy) - currentU.getData(i,j,myNx,myNy) * currentU.getData(i,j,myNx,myNy) * currentV.getData(i,j,myNx,myNy));
            nextV.updateData(i, j, myNx, myNy, updateV); 
         }
    }

    for (int i = 1; i <= myNy; i++)
    {
        for (int j = 1; j <= myNx; j++)
        {
            double newU = nextU.getData(i,j,myNx,myNy);
            double newV = nextV.getData(i,j,myNx,myNy);


            currentU.updateData(i, j, myNx, myNy, newU); 
            currentV.updateData(i, j, myNx, myNy, newV); 
        }
    }

    currentU.haloExchange(NumColumns, N);
    currentV.haloExchange(NumColumns, N);

}

State Brusselator::evaluateRHS(double coeff)
{

    ParMesh& currentU = myVars.getField(0).getCurrent();
    ParMesh& nextU    = myVars.getField(0).getNext();

    ParMesh& currentV = myVars.getField(1).getCurrent();
    ParMesh& nextV    = myVars.getField(1).getNext();

    int NumColumns = myVars.getField(0).getCols();
    int N          = myVars.getField(0).getN();

    int numRows = currentU.getSize() / NumColumns;
    int myNx = N / NumColumns;
    int myNy = N / numRows;
    int stride = myNx + 2;     // row width including left/right halos

    double dx = myVars.getField(0).getDomain() / N;
    double dy = myVars.getField(0).getDomain() / N;



    // ---------------- INTERIOR ROWS WITH SIDE HALOS ----------------
    for (int i = 1; i <= myNy; i++)
    {
        for (int j = 1; j <= myNx; j++)
        {
            double Diff_U = myParams.Du;
            double Diff_V = myParams.Dv;
            double A = myParams.A;
            double B = myParams.B;

            double UnewDiff = Diff_U * ( (currentU.getData(i+1,j,myNx,myNy) - 2 * currentU.getData(i,j,myNx,myNy) + currentU.getData(i-1,j,myNx,myNy)) / (dx*dx) + (currentU.getData(i,j+1,myNx,myNy) - 2 * currentU.getData(i,j,myNx,myNy) + currentU.getData(i,j-1,myNx,myNy)) / (dy*dy));
            double VnewDiff = Diff_V * ( (currentV.getData(i+1,j,myNx,myNy) - 2 * currentV.getData(i,j,myNx,myNy) + currentV.getData(i-1,j,myNx,myNy)) / (dx*dx) + (currentV.getData(i,j+1,myNx,myNy) - 2 * currentV.getData(i,j,myNx,myNy) + currentV.getData(i,j-1,myNx,myNy)) / (dy*dy));

            double updateU =  coeff * (UnewDiff + A + currentU.getData(i,j,myNx,myNy) * currentU.getData(i,j,myNx,myNy) * currentV.getData(i,j,myNx,myNy) - (B+1)*currentU.getData(i,j,myNx,myNy));
            nextU.updateData(i, j, myNx, myNy, updateU); 
        
            double updateV =  coeff * (VnewDiff + B*currentU.getData(i,j,myNx,myNy) - currentU.getData(i,j,myNx,myNy) * currentU.getData(i,j,myNx,myNy) * currentV.getData(i,j,myNx,myNy));
            nextV.updateData(i, j, myNx, myNy, updateV); 
         }
    }

    for (int i = 1; i <= myNy; i++)
    {
        for (int j = 1; j <= myNx; j++)
        {
            double newU = nextU.getData(i,j,myNx,myNy);
            double newV = nextV.getData(i,j,myNx,myNy);


            currentU.updateData(i, j, myNx, myNy, newU); 
            currentV.updateData(i, j, myNx, myNy, newV); 
        }
    }

    currentU.haloExchange(NumColumns, N);
    currentV.haloExchange(NumColumns, N);

    return myVars; 

}

void Brusselator::initialize(std::mt19937 gen, std::uniform_real_distribution<> Distribution)
{
    ParMesh& currentU = myVars.getField(0).getCurrent();
    ParMesh& nextU    = myVars.getField(0).getNext();

    ParMesh& currentV = myVars.getField(1).getCurrent();
    ParMesh& nextV    = myVars.getField(1).getNext();

    int NumColumns = myVars.getField(0).getCols();
    int N          = myVars.getField(0).getN();

    int numRows = currentU.getSize() / NumColumns;
    int myNx = N / NumColumns;
    int myNy = N / numRows;
    int stride = myNx + 2;     // row width including left/right halos

    double dx = myVars.getField(0).getDomain() / N;
    double dy = myVars.getField(0).getDomain() / N;

    for (int i = 1; i <= myNy; i++)
    {
        for (int j = 1; j <= myNx; j++)
        {

            double updateU = Distribution(gen) + myParams.A;
            double updateV = Distribution(gen) + myParams.B;

            currentU.updateData(i, j, myNx, myNy, updateU);
            currentV.updateData(i, j, myNx, myNy, updateV); 
        }
    }


    currentU.haloExchange(NumColumns, N);
    currentV.haloExchange(NumColumns, N);

}

void Brusselator::setState(State& newState)
{
    int numFields = myVars.getNumdFields();

    for (int k = 0; k < numFields; k++)
    {
        // Source data
        const std::vector<double>& src =
            newState.getField(k).getCurrent().getRawData();

        // Destination data
        std::vector<double>& dest =
            myVars.getField(k).getCurrent().getRawData();

        int size = src.size();

        for (int i = 0; i < size; i++)
        {
            dest[i] = src[i];
        }
    }
}
