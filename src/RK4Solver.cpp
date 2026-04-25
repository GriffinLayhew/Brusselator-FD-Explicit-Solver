#include "RK4Solver.hpp"

// Constructor Definition
RK4Solver::RK4Solver(Brusselator& Problem, double dt, int Steps, int SF) : myProblem(Problem), timeStep(dt), numSteps(Steps), saveFreq(SF), myOutput(myProblem.getVars(), "TimeStep.out") 
{

}
// Non-Modifying Definitions


// Modifying Member Functions
void RK4Solver::doNextStep()
{
    double alpha1 = timeStep / 4.0;
    double alpha2 = timeStep / 3.0;
    double alpha3 = timeStep / 2.0;
    double alpha4 = timeStep / 1.0;

    State CurrentState   = myProblem.getVars();
    State OriginalState  = myProblem.getVars();

    CurrentState         = OriginalState + myProblem.evaluateRHS(alpha1);
    myProblem.setState                                    (CurrentState);
    CurrentState         = OriginalState + myProblem.evaluateRHS(alpha2);
    myProblem.setState                                    (CurrentState);
    CurrentState         = OriginalState + myProblem.evaluateRHS(alpha3);
    myProblem.setState                                    (CurrentState);
    CurrentState         = OriginalState + myProblem.evaluateRHS(alpha4);
    myProblem.setState                                    (CurrentState);
}

void RK4Solver::Solve()
{
    for (int i=0; i<numSteps; i++)
    {
        this->doNextStep();

        if (i % saveFreq == 0)
        {
            myOutput.printToFile();
        }
    }
}



// Helper Definitions

