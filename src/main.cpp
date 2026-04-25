/*
    AUTHOR:                     Griffin Layhew
    DATE:                       04/22/26
    PROGRAM TITLE:              Brusselator Reaction
    PROGRAM DESCRIPTION:
            " This program... "


*/


#include "Utilities.hpp"
#include "ParMesh.hpp"
#include "Brusselator.hpp"
#include "RK4Solver.hpp"


int main(int argc, char* argv[])
{
    /*
    =============================================
                    Intiialize MPI
    =============================================
    */
    int rank, size;
    MPI_Status status; 

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    /*
    =============================================
                Intiialize Constants
    =============================================
    */

    int N = 128;

    double Domain = 1.0;
    int NumColumns = 2;

    Parameters params;
    params.A = 1.0;
    params.B = 3.0;
    params.Du = 5.0e-5;
    params.Dv = 5.0e-6;

    double      timeStep = 0.1;
    int         numSteps = 100000;
    int         saveFreq = 1000;

    /*
    =============================================
    Intiialize Random Number Generator with Set Seed
    =============================================
    */

    std::uniform_real_distribution<> Distribution(-1,1);

    /*
    =============================================
              Initialize Field Data
    =============================================
    */


    // Testing

    std::vector<std::string> varNames = {"U", "V"};
    State theState(rank, size, N, Domain, NumColumns, varNames);
    Brusselator Problem(params, theState);
    Problem.initialize(Distribution, rank);

    RK4Solver Solver(Problem, timeStep, numSteps, saveFreq);
    Solver.Solve();


    if (rank == 3)
    {
      Problem.display();
    }


    MPI_Finalize();
    return 0;
}