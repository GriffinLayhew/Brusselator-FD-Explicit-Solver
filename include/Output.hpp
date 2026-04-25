#pragma once

#include "Utilities.hpp"
#include "State.hpp"

class Output
{
    public:
    // Constructor
        Output(State& aState, std::string aFileName) : myState(aState), fileName(aFileName) {};



    // Non-Modifying Member Functions
        void printToFile();


    // Modifying Member Functions



    private:
        State& myState;
        std::string  fileName;




};

// Helper Functions

