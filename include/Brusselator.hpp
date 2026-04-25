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
        State   evaluateRHS(double);
        void    forwardEuler(double);
        void    initialize(std::mt19937 gen, std::uniform_real_distribution<> Distribution);
        State&  getVars() { return myVars; }
        const   State& getVars() const { return myVars; }
        void    setState(State&);

    private:
        Parameters myParams;
        State&     myVars;

};

// Helper Functions


