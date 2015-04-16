//Euler made by Marnix van Schrojenstein Lantman

#include "Euler.h"

#include "Logger.h"

auto& dimension = Base::register_argument<std::size_t>('D', "dim", "number of dimensions in the problem");
auto& numOfElements = Base::register_argument<std::size_t>('n', "numElems", "number of elements per dimension", true);
auto& polynomialOrder = Base::register_argument<std::size_t>('p', "order", "polynomial order of the solution", true);

auto& numOfOutputFrames = Base::register_argument<std::size_t>('O', "numOfOutputFrames", "Number of frames to output", false, 200);
auto& startTime = Base::register_argument<double>('S', "startTime", "start time of the simulation", false, 0.0);
auto& endTime = Base::register_argument<double>('T', "endTime", "end time of the simulation", false, 5.0);
auto& dt = Base::register_argument<double>('d', "timeStepSize", "time step of the simulation", false, 0.001);

int main (int argc, char **argv){

	Base::parse_options(argc, argv);

	logger(WARN,"WARNING: Timestep is determined a priori. Stability Criteria might not be satisfied!");
    // Set parameters for the PDE.
    const Base::MeshType meshType = Base::MeshType::RECTANGULAR;
    const Base::ButcherTableau * const ptrButcherTableau = Base::AllTimeIntegrators::Instance().getRule(3,3,true);

    //Set variable names and number of parameters
    std::vector<std::string> variableNames;
    variableNames.push_back("q1");

    for(std::size_t i = 1; i <= dimension.getValue()+1; i++) // +1 is for the energy equation
    {
        std::string variableName = "q" + std::to_string(i+1);
        variableNames.push_back(variableName);
    }

   // Calculate number of variables
   const std::size_t numOfVariables = 2+dimension.getValue();

    // Create problem solver 'test', that can solve the euler equations.
    Euler test(dimension.getValue(), numOfVariables, polynomialOrder.getValue(), ptrButcherTableau);

    // Create the mesh, a simple square domain
    test.createMesh(numOfElements.getValue(), meshType);

    // Solve the problem over time interval [startTime,endTime].
    test.solve(startTime.getValue(), endTime.getValue(), dt.getValue(), numOfOutputFrames.getValue(), true);

    return 0;
}
