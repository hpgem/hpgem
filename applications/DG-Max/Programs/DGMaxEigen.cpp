
#include <exception>
#include "Base/CommandLineOptions.h"

#include "DGMaxLogger.h"
#include "DGMaxProgramUtils.h"
#include "Algorithms/DivDGMaxEigenValue.h"


// File name of the mesh file, e.g. -m mesh.hpgem
auto& meshFile = Base::register_argument<std::string>('m', "meshFile",
        "The hpgem meshfile to use", true);
// Polynomial order of the basis functions, e.g. -p 2
auto& p = Base::register_argument<std::size_t>('p', "order",
        "Polynomial order of the solution", true);
// Number of eigenvalues to compute, e.g. -e 40
auto& numEigenvalues = Base::register_argument<std::size_t>('e', "eigenvalues",
        "The number of eigenvalues to compute", false, 24);

// Compute a single point --point 1,0.5,0
// Which corresponds to the point pi, 0.5pi, 0 in k-space
auto& pointMode = Base::register_argument<std::string>('\0', "point",
        "Compute only this single k-space-point", false);

// Number of steps to take in the k-space walk along the cube edge
// e.g. --steps 20
auto& steps = Base::register_argument<std::size_t>('\0', "steps",
        "Steps for the k-space walk", false, 10);

// Penalty parameter, e.g. -a b5b0b5
auto& pparams = Base::register_argument<std::string>('a', "penalty",
        "Penalty parameters and fluxes to use, e.g. b5b0b5.0\n"
        "\t  The letters are for the flux types and can be either b (Brezzi) "
        "or i (Interior penalty). The numbers are for the three penalty values.", false);

// Dimension, e.g. -d 2
auto& d = Base::register_argument<std::size_t>('d', "dimension",
        "The dimension of the problem", true);

template<std::size_t DIM>
void runWithDimension();
template<std::size_t DIM>
typename DivDGMaxDiscretization<DIM>::Stab parsePenaltyParmaters();
template<std::size_t DIM>
KSpacePath<DIM> parsePath();

int main(int argc, char** argv)
{
    Base::parse_options(argc, argv);
    initDGMaxLogging();
    DGMax::printArguments(argc, argv);

    time_t start, end;
    time(&start);

    const std::size_t dimension = d.getValue();
    try
    {
        switch(dimension)
        {
            case 2:
                runWithDimension<2>();
                break;
            case 3:
                runWithDimension<3>();
                break;
            default:
                logger.assert_always(false, "Invalid dimension %", dimension);
        }
    }
    catch (const std::exception& e)
    {
        DGMaxLogger(ERROR, e.what());
        exit(2);
    }
    catch (const char* message)
    {
        DGMaxLogger(ERROR, message);
        exit(1);
    }
    time(&end);
    DGMaxLogger(INFO, "Runtime %s", end - start);
    return 0;
}

template<std::size_t DIM>
void runWithDimension()
{
    typename DivDGMaxDiscretization<DIM>::Stab stab = parsePenaltyParmaters<DIM>();

    // 2 unknowns, 1 time level
    Base::ConfigurationData configData(2, 1);
    auto mesh = DGMax::readMesh<DIM>(meshFile.getValue(), &configData);
    DivDGMaxEigenValue<DIM> solver(*mesh);
    // TODO: Parameterize
    KSpacePath<DIM> path = parsePath<DIM>();
    EigenValueProblem<DIM> input(path, numEigenvalues.getValue());
    typename DivDGMaxEigenValue<DIM>::Result result = solver.solve(input, stab, p.getValue());
    if (Base::MPIContainer::Instance().getProcessorID() == 0)
    {
        result.printFrequencies();
        result.writeFrequencies("frequencies.csv");
    }
}

template<std::size_t DIM>
KSpacePath<DIM> parsePath()
{
    if(pointMode.isUsed())
    {
        const std::string& pointString = pointMode.getValue();
        typename KSpacePath<DIM>::KPoint point;
        std::size_t start = 0;
        // Parse DIM comma separated values
        for(std::size_t i = 0; i < DIM; ++i)
        {
            if(start >= pointString.size())
            {
                throw std::invalid_argument("Not enough coordinates for a reciprocal point");
            }
            std::size_t end = pointString.find_first_of(',', start);
            if(end == std::string::npos)
            {
                end = pointString.length();
            }
            point[i] = std::stod(pointString.substr(start, end - start));
            start = end + 1;// Skip the comma itself
        }
        if(start != pointString.length() + 1)
        {
            throw std::invalid_argument("Too many coordinates for a reciprocal point");
        }
        // Treat it as reduced coordinate.
        point *= M_PI;
        return KSpacePath<DIM>::singleStepPath(point);
    }
    else
    {
        if(!steps.isUsed())
        {
            logger(INFO, "Using default number of steps %", steps.getValue());
        }
        return KSpacePath<DIM>::cubePath(steps.getValue(), false);
    }
}

template<std::size_t DIM>
typename DivDGMaxDiscretization<DIM>::Stab parsePenaltyParmaters()
{
    if(pparams.isUsed())
    {
        typename DivDGMaxDiscretization<DIM>::Stab stab;
        std::string input = pparams.getValue();
        std::vector<bool> useBrezzi;
        std::vector<double> values;

        std::size_t index = 0;
        bool error = false;
        // Only need three parameters
        while(input.size() > index && !error && useBrezzi.size() < 3)
        {
            // Parse the flux type
            char fluxType = input[index++];
            if(fluxType == 'b')
                useBrezzi.emplace_back(true);
            else if(fluxType == 'i')
                useBrezzi.emplace_back(false);
            else
            {
                DGMaxLogger(ERROR, "Unknown flux type %", fluxType);
                error = true;
                break;
            }
            // Find the numeric digits
            std::size_t startIndex = index;
            while(input.size() > index)
            {
                char c = input[index];
                if(std::isdigit(c) || c == '.')
                {
                    index++;
                }
                else
                {
                    // Not a valid part of the number
                    break;
                }
            }
            // Parse the number (if present)
            if(index != startIndex)
            {
                // There is some numeric content
                double value = std::stod(input.substr(startIndex, index - startIndex));
                values.emplace_back(value);
            }
            else if(index == input.size())
            {
                DGMaxLogger(ERROR, "Not enough input for the stabilization parameter");
                error = true;
            }
            else
            {
                DGMaxLogger(ERROR, "No number at position % of the stabilization parameters", index);
                error = true;
            }
        }

        // Check the validity of the result
        if(!error && useBrezzi.size() < 3)
        {
            DGMaxLogger(ERROR, "Not enough stabilization parameters only parsed %", useBrezzi.size());
            error = true;
        }
        if(!error && index != input.size())
        {
            DGMaxLogger(ERROR, "Unconsumed stabilization parameter input at %", index);
            error = true;
        }
        if(!error)
        {
            typename DivDGMaxDiscretization<DIM>::Stab result;
            result.stab1 = values[0];
            result.stab2 = values[1];
            result.stab3 = values[2];

            using FLUX = typename DivDGMaxDiscretization<DIM>::FluxType;

            result.fluxType1 = useBrezzi[0] ? FLUX::BREZZI : FLUX::IP;
            result.fluxType2 = useBrezzi[1] ? FLUX::BREZZI : FLUX::IP;
            result.fluxType3 = useBrezzi[2] ? FLUX::BREZZI : FLUX::IP;
            return result;
        }
        else
        {
            throw std::invalid_argument("Invalid stabilization parameter");
        }
    }
    else
    {
        // Default values
        typename DivDGMaxDiscretization<DIM>::Stab stab;
        stab.stab1 = 5;
        stab.stab2 = 0;
        stab.stab3 = 5;
        stab.setAllFluxeTypes(DivDGMaxDiscretization<DIM>::FluxType::BREZZI);
        return stab;
    }
}