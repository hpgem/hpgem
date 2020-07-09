
#include "BandstructureUtils.h"

#include "LinearAlgebra/SmallMatrix.h"

using namespace hpgem;


template <std::size_t DIM>
bool next(std::array<int, DIM> &arr, int min, int max) {
    for (std::size_t i = 0; i < DIM; ++i) {
        if (arr[i] < max) {
            arr[i]++;
            return true;
        }
        arr[i] = min;
    }
    return false;
}

double intervalDist(double x, double xmin, double xmax) {
    logger.assert_debug(xmin <= xmax, "xmin should not be more than xmax");
    if (x < xmin) {
        return xmin - x;
    }
    if (x < xmax) {
        return 0;
    } else {
        return x - xmax;
    }
}

std::map<double, std::size_t> group(std::vector<double> vect,
                                    double tolerance) {
    std::map<double, std::size_t> result;
    {
        auto iter = vect.begin();
        while (iter != vect.end()) {
            double entry = *iter;
            std::size_t count = 1;
            iter++;
            while (iter != vect.end() && std::abs(*iter - entry) < tolerance) {
                count++;
                iter++;
            }
            result[entry] = count;
        }
    }
    return result;
}

template <std::size_t DIM>
std::vector<LatticePoint<DIM>> LatticePoint<DIM>::getNeighbours() const {
    std::array<int, DIM> offset, ncoords;
    offset.fill(-1);

    std::vector<LatticePoint> result;
    result.reserve(std::round(std::pow(3, DIM) - 1));

    do {
        bool zeroOffset = true;
        for (std::size_t i = 0; i < DIM; ++i) {
            zeroOffset &= offset[i] == 0;
            ncoords[i] = coords_[i] + offset[i];
        }
        if (zeroOffset) {
            continue;
        }
        result.emplace_back(ncoords);
    } while (next(offset, -1, 1));
    return result;
}

template <std::size_t DIM>
LinearAlgebra::SmallVector<DIM> LatticePoint<DIM>::k(
    const LatticePoint::Basis &basis) const {
    LinearAlgebra::SmallVector<DIM> result;
    result.set(0);
    for (std::size_t i = 0; i < DIM; ++i) {
        result += coords_[i] * basis[i];
    }
    return result;
}

template <std::size_t DIM>
std::vector<LatticePoint<DIM>> LatticePoint<DIM>::getNearestNeighbours(
    const LinearAlgebra::SmallVector<DIM> &kpoint,
    const LatticePoint::Basis &reciprocalBasis) {
    // Insert initial points.
    std::vector<LatticePoint> result;
    result.reserve(1 << (DIM + 1));
    LinearAlgebra::SmallMatrix<DIM, DIM> basis;
    for (std::size_t vecId = 0; vecId < DIM; ++vecId) {
        for (std::size_t coordId = 0; coordId < DIM; ++coordId) {
            basis(coordId, vecId) = reciprocalBasis[vecId][coordId];
        }
    }
    LinearAlgebra::SmallVector<DIM> realCoords = kpoint;
    basis.solve(realCoords);
    // We now have the coordinates in terms of the reciprocal lattice.
    // So add the nearest neighbours to this point
    std::array<int, DIM> roundedCoords, offset, icoords;
    offset.fill(-1);
    for (std::size_t i = 0; i < DIM; ++i) {
        roundedCoords[i] = std::round(realCoords[i]);
    }
    do {
        for (std::size_t i = 0; i < DIM; ++i) {
            icoords[i] = roundedCoords[i] + offset[i];
        }
        result.emplace_back(icoords);
    } while (next(offset, -1, 1));
    return result;
}

template class LatticePoint<1>;
template class LatticePoint<2>;
template class LatticePoint<3>;