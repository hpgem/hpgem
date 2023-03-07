
#ifndef HPGEM_APP_BANDSTRUCTUREUTILS_H
#define HPGEM_APP_BANDSTRUCTUREUTILS_H

#include "LinearAlgebra/SmallVector.h"

#include <map>
#include <vector>

using namespace hpgem;

// Several utility functions for computing band structures

/// \brief Lexiographical next array.
///
/// \tparam DIM The dimension in the array
/// \param arr The array to modify to the next item.
/// \param min The minimum value of the entries
/// \param max The maximum value of the entries
/// \return Whether a next array is available
template <std::size_t DIM>
bool next(std::array<int, DIM>& arr, int min, int max);

/// \brief Compute distance to an interval
///
/// \param x The point
/// \param xmin The minumum coordinate of the interval
/// \param xmax The maximum coordinate of the interval
/// \return The distance between the point and the interval.
double intervalDist(double x, double xmin, double xmax);

/// \brief Group entries of a vector if close enough
///
/// \param vect The vector with entries
/// \param tolerance The tolerance below which elements should be grouped
/// \return A map from keys to counts.
std::map<double, std::size_t> group(std::vector<double> vect, double tolerance);

/// \brief Lattice coordinates
template <std::size_t DIM>
struct LatticePoint {
    using Basis = std::array<LinearAlgebra::SmallVector<DIM>, DIM>;

    LatticePoint() {
        for (std::size_t i = 0; i < DIM; ++i) {
            coords_[i] = 0;
        }
    }

    LatticePoint(std::array<int, DIM> coords) : coords_(coords) {}

    LatticePoint(const LatticePoint<DIM>& other) = default;
    LatticePoint(LatticePoint<DIM>&& other) noexcept = default;

    LatticePoint<DIM>& operator=(const LatticePoint<DIM>& other) = default;
    LatticePoint<DIM>& operator=(LatticePoint<DIM>&& other) noexcept = default;

    bool operator==(const LatticePoint<DIM>& other) const {
        return coords_ == other.coords_;
    }

    // Lexicographic order
    bool operator<(const LatticePoint<DIM>& other) const {
        return coords_ < other.coords_;
    }

    int operator[](int index) const { return coords_[index]; }

    std::vector<LatticePoint> getNeighbours() const;

    LinearAlgebra::SmallVector<DIM> k(const Basis& basis) const;

    /// \brief Get the lattice points close to a kpoint for a basis
    ///
    /// For stability not only the closest point is returned, but also its
    /// neighbours.
    /// \param kpoint The point in k space
    /// \param reciprocalBasis The basis of the reciprocal lattice
    /// \return The points in kspace that are close
    static std::vector<LatticePoint> getNearestNeighbours(
        const LinearAlgebra::SmallVector<DIM>& kpoint,
        const Basis& reciprocalBasis);

   private:
    std::array<int, DIM> coords_;
};

#endif  // HPGEM_APP_BANDSTRUCTUREUTILS_H
