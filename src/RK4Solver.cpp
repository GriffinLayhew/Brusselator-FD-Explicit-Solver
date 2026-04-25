#include "RK4Solver.hpp"

// Constructor Definition
RK4Solver::RK4Solver(Brusselator& Problem, double dt, int Steps, int SF) : myProblem(Problem), timeStep(dt), numSteps(Steps), saveFreq(SF), myOutput(myProblem.getVars(), "TimeStep.out") 
{

}
// Non-Modifying Definitions


// Modifying Member Functions
void RK4Solver::doNextStep()
{
    State Fn = myProblem.getVars();

    State F1    = Fn + myProblem.evaluateRHS(Fn) * (timeStep / 4.0);
    State F2    = Fn + myProblem.evaluateRHS(F1) * (timeStep / 3.0);
    State F3    = Fn + myProblem.evaluateRHS(F2) * (timeStep / 2.0);
    State Fnew  = Fn + myProblem.evaluateRHS(F3) * timeStep;

    myProblem.setState(Fnew);
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

