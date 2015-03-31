//Euler made by Marnix van Schrojenstein Lantman

#include "Euler.h"

#include "Logger.h"

int main (int argc, char **argv){

	Base::parse_options(argc, argv);

	logger(WARN,"WARNING: Timestep is determined a priori. CLF condition might not be satisfied!");
    logger(WARN,"WARNING: Pressure might be incorrectly calculated");
	logger(INFO,"INFO: test.");
    // Set parameters for the PDE.
    const std::size_t dimension = 2;
    const std::size_t numOfElements = 20;
    const std::size_t polynomialOrder = 0;
    const Base::MeshType meshType = Base::MeshType::TRIANGULAR;
    const Base::ButcherTableau * const ptrButcherTableau = Base::AllTimeIntegrators::Instance().getRule(1,1);

    //Set variable names and number of parameters
    std::vector<std::string> variableNames;
    variableNames.push_back("q1");

    for(std::size_t i = 1; i <= dimension+1; i++) // +1 is for the energy equation
    {
        std::string variableName = "q" + std::to_string(i+1);
        variableNames.push_back(variableName);
    }

   // Calculate number of variables
   const std::size_t numOfVariables = 2+dimension;

    // Create problem solver 'test', that can solve the euler equations.
    Euler test(dimension, numOfVariables, polynomialOrder, ptrButcherTableau);

    // Create the mesh, a simple square domain
    test.createMesh(numOfElements, meshType);

    const double startTime = 0.0;
    const double endTime = 5.0;
    double dt = 0.0001;
    const std::size_t numOfOutputFrames = 100;

    // Solve the problem over time interval [startTime,endTime].
    test.solve(startTime, endTime, dt, numOfOutputFrames, false);

    return 0;
}
