#pragma once

#include "Utilities.hpp"
#include "Brusselator.hpp"
#include "Output.hpp"

class RK4Solver
{
    public:
    // Constructor
        RK4Solver(Brusselator& Problem, double dt, int Steps, int SF);


    // Non-Modifying Member Functions


    // Modifying Member Functions
        void doNextStep();
        void Solve();



    private:
        Brusselator& myProblem;
        double       timeStep;
        int          numSteps;
        int          saveFreq;
        Output       myOutput;



};

// Helper Functions

