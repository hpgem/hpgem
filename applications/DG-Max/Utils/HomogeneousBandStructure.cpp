
#include "HomogeneousBandStructure.h"

#include "LinearAlgebra/SmallMatrix.h"

#include <set>

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
struct LatticePoint
{
    double frequency_;
    std::array<int, DIM> coords_;

    LatticePoint(std::array<int, DIM> coords,
            const std::array<LinearAlgebra::SmallVector<DIM>, DIM>& lattice,
            const LinearAlgebra::SmallVector<DIM>& kpoint, double permitivity)
            : coords_ (coords)
    {
        LinearAlgebra::SmallVector<DIM> latticePoint;
        for(std::size_t i = 0; i < DIM; ++i)
        {
            latticePoint += lattice[i] * coords[i];
        }
        // omega = c |klat - kpoint| / n = \klat - kpoint| / sqrt(mu eps eps_r)
        // where we take eps = 1, mu = 1
        frequency_ = (kpoint - latticePoint).l2Norm() / std::sqrt(permitivity);
    }

    LatticePoint(const LatticePoint<DIM>& other) = default;
    LatticePoint(LatticePoint<DIM>&& other) noexcept = default;

    void addNeighbours(std::set<LatticePoint<DIM>>& queue,
            const std::array<LinearAlgebra::SmallVector<DIM>, DIM>& lattice,
            const LinearAlgebra::SmallVector<DIM>& kpoint,
            double permittivity)
    {
        std::array<int, DIM> offsets;
        std::array<int, DIM> newCoords;
        offsets.fill(-1);
        // Iterate over all offsets.
        do
        {
            // We need an offset, not offsets = {0,..,0}, which is the same as
            // the current lattice point
            bool allZero = true;
            for(std::size_t i = 0; i < DIM; ++i)
            {
                allZero &= offsets[i] == 0;
                newCoords[i] = offsets[i] + coords_[i];
            }
            if(allZero)
            {
                continue;
            }
            LatticePoint<DIM> toAdd(newCoords, lattice, kpoint, permittivity);
            if(toAdd.frequency_ > frequency_)
            {
                // Prevent adding the same point multiple times
                // Could be improved by using frequency_
                auto i = queue.begin();
                for(; i != queue.end(); ++i)
                {
                    if(i->matchesCoordinates(toAdd))
                        break;
                }
                if(i == queue.end())
                {
                    queue.insert(toAdd);
                }
            }
        } while(next(offsets, -1, 1));
    }

    bool matchesCoordinates(LatticePoint<DIM>& other) const
    {
        for(std::size_t i = 0; i < DIM; ++i)
        {
            if(coords_[i] != other.coords_[i])
                return false;
        }
        return true;
    }

    bool operator < (const LatticePoint<DIM>& other) const
    {
        if (frequency_ != other.frequency_)
        {
            return frequency_ < other.frequency_;
        }
        // Coordinate wise
        for(std::size_t i = 0; i < DIM; ++i)
        {
            if(coords_[i] != other.coords_[i])
                return coords_[i] < other.coords_[i];
        }
        return false;
    }
};

template<std::size_t DIM>
std::map<double, std::size_t> HomogeneousBandStructure<DIM>::computeSpectrum(
        LinearAlgebra::SmallVector<DIM> kpoint, int numberOfModes)
{
    std::set<LatticePoint<DIM>> queue;
    {
        // Insert initial points.
        LinearAlgebra::SmallMatrix<DIM, DIM> basis;
        for(std::size_t vecId = 0; vecId < DIM; ++vecId)
        {
            for(std::size_t coordId = 0; coordId < DIM; ++coordId)
            {
                basis(coordId, vecId) = reciprocalVectors_[vecId][coordId];
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
            queue.emplace(icoords, reciprocalVectors_, kpoint, permittivity_);
        } while(next(offset, -1, 1));
    }


    // Coordinates of the reciprocal lattice point with respect to the
    // reciprocal lattice basis of reciprocalVectors_.
    std::vector<double> frequencies;
    frequencies.reserve(numberOfModes + 16ul);

    double maxFrequency = 0;
    while(!queue.empty())
    {
        // Retrieve the lowest element
        auto next = queue.begin();
        LatticePoint<DIM> point = *next;
        queue.erase(next);

        // To ensure we get the right multiplicity we need to keep going even if
        // we already added 'numberOfModes'. However, as queue is ordered by
        // frequency, we can stop if we go beyond maxFrequency.
        if(frequencies.size() >= numberOfModes && point.frequency_ > maxFrequency + 1e-5)
        {
            // Everything added.
            break;
        }
        // Insert in the frequency table
        double frequency = point.frequency_;
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
            point.addNeighbours(queue, reciprocalVectors_, kpoint, permittivity_);
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

template class HomogeneousBandStructure<2>;
template class HomogeneousBandStructure<3>;

