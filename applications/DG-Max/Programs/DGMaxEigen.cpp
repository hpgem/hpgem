
#include <exception>
#include "Base/CommandLineOptions.h"

#include "DGMaxLogger.h"
#include "DGMaxProgramUtils.h"
#include "Algorithms/DivDGMaxEigenValue.h"
#include "Algorithms/DGMaxEigenValue.h"


// File name of the mesh file, e.g. -m mesh.hpgem
auto& meshFile = Base::register_argument<std::string>('m', "meshFile",
        "The hpgem meshfile to use", true);
// Polynomial order of the basis functions, e.g. -p 2
auto& p = Base::register_argument<std::size_t>('p', "order",
        "Polynomial order of the solution", true);
// Number of eigenvalues to compute, e.g. -e 40
auto& numEigenvalues = Base::register_argument<std::size_t>('e', "eigenvalues",
        "The number of eigenvalues to compute", false, 24);

auto& method = Base::register_argument<std::string>('\0', "method",
        "The method to be used, either 'DGMAX', 'DGMAX_PROJECT' or 'DIVDGMAX' (default)",
        false, "DIVDGMAX");

// Compute a single point --point 1,0.5,0 or a path of points
// [steps@]0,0:1,0:1,1
// Which corresponds to the point pi, 0.5pi, 0 in k-space and a path from 0,0
// via pi,0 to pi,pi in kspace taking each time "steps" steps, for each line,
// (excluding first point, including last).
// Untested with anything else than 0.0 as the first point of a path
auto& pointMode = Base::register_argument<std::string>('\0', "points",
        "Compute a single point in or a set of lines through k-space", false);

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

auto& structure = Base::register_argument<std::size_t>('\0', "structure",
        "Structure to use", false, 0);

template<std::size_t DIM>
void runWithDimension();
double parseDGMaxPenaltyParameter();
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
    bool useDivDGMax = true;
    bool useProjector = false;
    std::size_t numberOfElementMatrices = 2;
    std::size_t unknowns = 0;
    // 2 unknowns, 1 time level
    if (method.getValue() == "DGMAX")
    {
        useDivDGMax = false;
        unknowns = 1;
    }
    else if (method.getValue() == "DGMAX_PROJECT")
    {
        useDivDGMax = false;
        useProjector = true;
        unknowns = 2;
        // 1 more is needed for the projector operator
        numberOfElementMatrices = 3;
    }
    else if (method.getValue() == "DIVDGMAX")
    {
        useDivDGMax = true;
        unknowns = 2;
    }
    else
    {
        logger(ERROR, "Invalid method {}, should be either DGMAX, DGMAX_PROJECT or DIVDGMAX", method.getValue());
        return;
    }

    Base::ConfigurationData configData(unknowns, 1);
    auto mesh = DGMax::readMesh<DIM>(meshFile.getValue(), &configData, [&](const Geometry::PointPhysical<DIM>& p) {
        // TODO: Hardcoded structure
        return jelmerStructure(p, structure.getValue());
    }, numberOfElementMatrices);
    logger(INFO, "Loaded mesh % with % local elements", meshFile.getValue(), mesh->getNumberOfElements());
    // TODO: Parameterize
    KSpacePath<DIM> path = parsePath<DIM>();
    EigenValueProblem<DIM> input(path, numEigenvalues.getValue());
    // Method dependent solving
    if (useDivDGMax)
    {
        DivDGMaxEigenValue<DIM> solver(*mesh);
        typename DivDGMaxDiscretization<DIM>::Stab stab = parsePenaltyParmaters<DIM>();
        typename DivDGMaxEigenValue<DIM>::Result result = solver.solve(input, stab, p.getValue());
        if (Base::MPIContainer::Instance().getProcessorID() == 0)
        {
            result.printFrequencies();
            result.writeFrequencies("frequencies.csv");
        }
    }
    else
    {
        DGMaxEigenValue<DIM> solver (*mesh, p.getValue(), useProjector);
        const double stab = parseDGMaxPenaltyParameter();
        DGMaxEigenvalueBase::SolverConfig config;
        config.stab_ = stab;
        config.useHermitian_ = true;
        typename DGMaxEigenValue<DIM>::Result result = solver.solve(input, config);
        if (Base::MPIContainer::Instance().getProcessorID() == 0)
        {
            result.printFrequencies();
            result.writeFrequencies("frequencies.csv");
        }
    }
}

/// Parse DIM comma separated numbers as the coordinates of a point.
/// \tparam DIM The dimension of the point
/// \param pointString The string containing the point coordinates
/// \param start The starting index in pointString
/// \param point The point (out)
/// \return The first index in pointString after the number.
template<std::size_t DIM>
std::size_t parsePoint(const std::string& pointString, std::size_t start, LinearAlgebra::SmallVector<DIM>& point)
{
    for(std::size_t i = 0; i < DIM; ++i)
    {
        if(start >= pointString.size())
        {
            throw std::invalid_argument("Not enough coordinates for a reciprocal point");
        }
        std::size_t len = 0;
        point[i] = std::stod(pointString.substr(start), &len);
        if(len == 0)
        {
            throw std::invalid_argument("No value parsed");
        }
        start += len;
        if(i < DIM - 1)
        {
            // Skip the comma, space, whatever that ended the point
            start++;
        }
    }
    return start;
}

template<std::size_t DIM>
KSpacePath<DIM> parsePath()
{
    if(pointMode.isUsed())
    {
        DGMax::PointPath<DIM> path = DGMax::parsePath<DIM>(pointMode.getValue());
        // Compensate for factor of pi in the reciprocal lattice
        for(std::size_t i = 0; i < path.points_.size(); ++i)
        {
            path.points_[i] *= M_PI;
        }
        // Default steps to 1.
        if(path.steps_ < 0)
        {
            path.steps_ = 1;
        }
        return KSpacePath<DIM>(path.points_, (std::size_t) path.steps_);
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

double parseDGMaxPenaltyParameter()
{
    if (pparams.isUsed())
    {
        std::size_t idx;
        try
        {
            double value = std::stod(pparams.getValue(), &idx);
            if (idx != pparams.getValue().size())
            {
                throw std::invalid_argument("Invalid stabilization parameter, should be a single number for DGMAX");
            }
            return value;
        }
        catch (const std::invalid_argument&)
        {
            throw std::invalid_argument("Invalid stabilization parameter, should be a single number for DGMAX");
        }
    }
    else
    {
        // Default
        return 100;
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
                if(std::isdigit(c) || c == '.' || c == '-' || c == 'e')
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
            logger(INFO, "Using fluxes and stabilization: %", result);
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