#pragma once

#include "Utilities.hpp"
#include "State.hpp"

typedef struct 
{  
    double A;
    double B;
    double Du;
    double Dv;
    
} Parameters;


class Brusselator
{
    public:
    // Constructor
        Brusselator(Parameters params, State& myState);

    // Non-Modifying Member Functions
        void    display();


    // Modifying Member Functions
        State   evaluateRHS(const State& inputState);
        void    forwardEuler(double);
        void    initialize(std::uniform_real_distribution<>& Distribution, int);
        State&  getVars() { return myVars; }
        const   State& getVars() const { return myVars; }
        void    setState(const State& newState);

    private:
        Parameters myParams;
        State&     myVars;

};

// Helper Functions


