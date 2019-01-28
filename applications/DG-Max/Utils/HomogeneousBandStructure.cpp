
#include "HomogeneousBandStructure.h"

#include "LinearAlgebra/SmallMatrix.h"

#include <queue>

template<std::size_t DIM>
HomogeneousBandStructure<DIM>::HomogeneousBandStructure(
        std::array<LinearAlgebra::SmallVector<DIM>, DIM> reciprocalVectors, double permittivity)
        : reciprocalVectors_ (reciprocalVectors)
        , permittivity_ (permittivity)
{}

template<std::size_t DIM>
bool next(std::array<int, DIM>& arr, int min, int max)
{
    for(std::size_t i = 0; i < DIM; ++i)
    {
        if(arr[i] < max)
        {
            arr[i]++;
            return true;
        }
        else
        {
            arr[i] = min;
        }
    }
    return false;
}

template<std::size_t DIM>
using Basis = std::array<LinearAlgebra::SmallVector<DIM>, DIM>;

template<std::size_t DIM>
struct LatticePoint
{
    LatticePoint() = default;

    LatticePoint(std::array<int, DIM> coords)
            : coords_ (coords)
    {}

    LatticePoint(const LatticePoint<DIM>& other) = default;
    LatticePoint(LatticePoint<DIM>&& other) noexcept = default;

    LatticePoint<DIM>& operator =(const LatticePoint<DIM>& other) = default;
    LatticePoint<DIM>& operator =(LatticePoint<DIM>&& other) noexcept = default;

    bool operator==  (const LatticePoint<DIM>& other) const
    {
        return coords_ == other.coords_;
    }

    // Lexicographic order
    bool operator < (const LatticePoint<DIM>& other) const
    {
        return coords_ < other.coords_;
    }

    std::vector<LatticePoint> getNeighbours() const
    {
        std::array<int, DIM> offset, ncoords;
        offset.fill(-1);

        std::vector<LatticePoint> result;
        result.reserve(std::round(std::pow(3, DIM) - 1));

        do
        {
            bool zeroOffset = true;
            for(std::size_t i = 0; i < DIM; ++i)
            {
                zeroOffset &= offset[i] == 0;
                ncoords[i] = coords_[i] + offset[i];
            }
            if (zeroOffset)
            {
                continue;
            }
            result.emplace_back(ncoords);
        }
        while(next(offset, -1, 1));
        return result;
    }

    LinearAlgebra::SmallVector<DIM> k(Basis<DIM>& basis) const
    {
        LinearAlgebra::SmallVector<DIM> result;
        result.set(0);
        for(std::size_t i = 0; i < DIM; ++i)
        {
            result += coords_[i] * basis[i];
        }
        return result;
    }

    static std::vector<LatticePoint> getNearestNeighbours(
            const LinearAlgebra::SmallVector<DIM>& kpoint,
            const Basis<DIM>& reciprocalBasis)
    {
        // Insert initial points.
        std::vector<LatticePoint> result;
        result.reserve(1 << (DIM + 1));
        LinearAlgebra::SmallMatrix<DIM, DIM> basis;
        for(std::size_t vecId = 0; vecId < DIM; ++vecId)
        {
            for(std::size_t coordId = 0; coordId < DIM; ++coordId)
            {
                basis(coordId, vecId) = reciprocalBasis[vecId][coordId];
            }
        }
        LinearAlgebra::SmallVector<DIM> coords = kpoint;
        basis.solve(coords);
        // We now have the coordinates in terms of the reciprocal lattice.
        // So add the nearest neighbours to this point
        std::array<int, DIM> rcoords, offset, icoords;
        offset.fill(-1);
        for(std::size_t i = 0; i < DIM; ++i)
        {
            rcoords[i] = std::round(coords[i]);
        }
        do {
            for(std::size_t i = 0; i < DIM; ++i)
            {
                icoords[i] = coords[i] + offset[i];
            }
            result.emplace_back(icoords);
        } while(next(offset, -1, 1));
        return result;
    }

private:
    std::array<int, DIM> coords_;
};


template<std::size_t DIM>
std::map<double, std::size_t> HomogeneousBandStructure<DIM>::computeSpectrum(
        LinearAlgebra::SmallVector<DIM> kpoint, int numberOfModes)
{
    // Using this order so that points are first ordered by frequency.
    using Key = std::tuple<double, LatticePoint<DIM>>;

    // Abusing a set as priority queue, so that we can check if points have
    // already been entered into it.
    std::set<Key> queue;
    {
        std::vector<LatticePoint<DIM>> points
                = LatticePoint<DIM>::getNearestNeighbours(kpoint, reciprocalVectors_);
        for(LatticePoint<DIM>& p : points)
        {
            LinearAlgebra::SmallVector<DIM> dk = p.k(reciprocalVectors_) - kpoint;
            double frequency = dk.l2Norm() / std::sqrt(permittivity_);
            queue.emplace(frequency, p);
        }
    }


    // Coordinates of the reciprocal lattice point with respect to the
    // reciprocal lattice basis of reciprocalVectors_.
    std::vector<double> frequencies;
    frequencies.reserve(numberOfModes + 16ul);

    double maxFrequency = 0;
    while(!queue.empty())
    {
        LatticePoint<DIM> point;
        double frequency;
        // Retrieve the lowest element
        auto next = queue.begin();
        std::tie(frequency, point) = *next;
        queue.erase(next);

        // To ensure we get the right multiplicity we need to keep going even if
        // we already added 'numberOfModes'. However, as queue is ordered by
        // frequency, we can stop if we go beyond maxFrequency.
        if(frequencies.size() >= numberOfModes && frequency > maxFrequency + 1e-5)
        {
            // Everything added.
            break;
        }
        // Insert in the frequency table
        maxFrequency = frequency;
        frequencies.emplace_back(frequency);
        // Two polarization modes for 3D
        if (DIM > 2)
        {
            frequencies.emplace_back(frequency);
        }
        // Zero mode has an extra mode, (thus for each cardinal direction once)
        if(std::abs(frequency) < 1e-12)
        {
            frequencies.emplace_back(frequency);
        }

        if(frequencies.size() < numberOfModes)
        {
            for(LatticePoint<DIM> n : point.getNeighbours())
            {
                double neighbourFrequency = (n.k(reciprocalVectors_) - kpoint).l2Norm()
                        / std::sqrt(permittivity_);
                if(neighbourFrequency > frequency)
                {
                    // Prevent double insertion by checking if the neighbour is
                    // already in the queue. We explicitly check the
                    // coordinates, as the frequency is a floating point number.
                    bool found = false;
                    for(auto iter = queue.begin(); iter != queue.end(); ++iter)
                    {
                        if(std::get<1>(*iter) == n)
                        {
                            found = true;
                            break;
                        }
                    }
                    if(!found)
                    {
                        queue.emplace(neighbourFrequency, n);
                    }
                }
            }
        }

    }
    // group the results
    std::map<double, std::size_t> result;
    {
        auto iter = frequencies.begin();
        while(iter != frequencies.end())
        {
            double frequency = *iter;
            std::size_t count = 1;
            iter++;
            while(iter != frequencies.end() && std::abs(*iter - frequency) < 1e-5)
            {
                count++;
                iter++;
            }
            result[frequency] = count;
        }
    }
    return result;
}

double intervalDist(double x, double xmin, double xmax)
{
    logger.assert_debug(xmin <= xmax, "xmin should not be more than xmax");
    if(x < xmin)
    {
        return xmin - x;
    }
    else if (x < xmax)
    {
        return 0;
    }
    else
    {
        return x - xmax;
    }
}

template<std::size_t DIM>
std::vector<typename HomogeneousBandStructure<DIM>::Line> HomogeneousBandStructure<DIM>::computeLines(
        LinearAlgebra::SmallVector<DIM> point1,
        LinearAlgebra::SmallVector<DIM> point2, double maxFrequency)
{
    const double l = (point1 - point2).l2Norm();
    LinearAlgebra::SmallVector<DIM> kdir = point2 - point1;
    kdir /= kdir.l2Norm();

    // Set of points that are considered (either in the queue or were in the
    // queue at some previous point).
    std::set<LatticePoint<DIM>> considered;
    // Queue of points that still have to be processed.
    std::queue<LatticePoint<DIM>> queue;
    // Identify the resulting lines/modes by the two coordinates, x,y
    // and associate with it the multiplicity.
    std::map<LinearAlgebra::SmallVector<2>, std::size_t> modes;

    for(LatticePoint<DIM> n : LatticePoint<DIM>::getNearestNeighbours(point1, reciprocalVectors_))
    {
        considered.emplace(n);
        queue.emplace(n);
    }
    while (!queue.empty())
    {
        LatticePoint<DIM> p = queue.front();
        queue.pop();
        // Compute the relevant parameters
        LinearAlgebra::SmallVector<DIM> dk = point1 - p.k(reciprocalVectors_);
        // See documentation in header for the interpretation of these values
        double x = dk * kdir;
        double y = (dk - x * kdir).l2Norm();
        // Compute the minimum frequency on this part of the line for filtering
        // lattice points which are above the maxFrequency
        // Note that x=0 corresponds to k1, and x=-l corresponds to k2.
        double xmin = intervalDist(x, -l, 0);
        double fmin = std::sqrt(xmin * xmin + y * y) / std::sqrt(permittivity_);
        if(fmin > maxFrequency)
        {
            continue;
        }
        // Add the point to the result
        LinearAlgebra::SmallVector<2> newMode ({x,y});
        bool added = false;
        for(auto& mode : modes)
        {
            if((mode.first - newMode).l2NormSquared() < 1e-10)
            {
                mode.second++;
                added = true;
            }
        }
        if(!added)
        {
            modes[newMode] = 1;
        }

        for(LatticePoint<DIM> n : p.getNeighbours())
        {
            // Enqueue those neighbours that have not already been visited.
            if(considered.find(n) != considered.end())
            {
                continue;
            }
            queue.push(n);
            considered.emplace(n);
        }
    }
    // construct the result
    std::vector<Line> result;
    for(auto& mode : modes)
    {
        result.emplace_back(l, mode.first[0], mode.first[1], permittivity_, mode.second);
    }
    return result;
}

template class HomogeneousBandStructure<2>;
template class HomogeneousBandStructure<3>;

