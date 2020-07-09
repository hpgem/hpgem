
#include "HomogeneousBandStructure.h"
#include "BandstructureUtils.h"

#include <queue>
#include <iostream>

using namespace hpgem;

template <std::size_t DIM>
HomogeneousBandStructure<DIM>::HomogeneousBandStructure(
    std::array<LinearAlgebra::SmallVector<DIM>, DIM> reciprocalVectors,
    double permittivity)
    : reciprocalVectors_(reciprocalVectors), permittivity_(permittivity) {}

template <std::size_t DIM>
std::map<double, std::size_t> HomogeneousBandStructure<DIM>::computeSpectrum(
    LinearAlgebra::SmallVector<DIM> kpoint, double maxFrequency) const {
    // Using this order so that points are first ordered by frequency.
    using Key = std::tuple<double, LatticePoint<DIM>>;

    // Abusing a set as priority queue, so that we can check if points have
    // already been entered into it.
    std::set<Key> queue;
    {
        std::vector<LatticePoint<DIM>> points =
            LatticePoint<DIM>::getNearestNeighbours(kpoint, reciprocalVectors_);
        for (LatticePoint<DIM>& p : points) {
            LinearAlgebra::SmallVector<DIM> dk =
                p.k(reciprocalVectors_) - kpoint;
            double frequency = dk.l2Norm() / std::sqrt(permittivity_);
            if (frequency > maxFrequency) continue;

            queue.emplace(frequency, p);
        }
    }

    // Coordinates of the reciprocal lattice point with respect to the
    // reciprocal lattice basis of reciprocalVectors_.
    std::vector<double> frequencies;

    while (!queue.empty()) {
        LatticePoint<DIM> point;
        double frequency;
        // Retrieve the lowest element
        auto next = queue.begin();
        std::tie(frequency, point) = *next;
        queue.erase(next);

        // Insert in the frequency table
        frequencies.emplace_back(frequency);
        // Two polarization modes for 3D
        if (DIM > 2) {
            frequencies.emplace_back(frequency);
        }
        // Zero mode has an extra mode, (thus for each cardinal direction once)
        if (std::abs(frequency) < 1e-12) {
            frequencies.emplace_back(frequency);
        }

        for (LatticePoint<DIM> n : point.getNeighbours()) {
            double neighbourFrequency =
                (n.k(reciprocalVectors_) - kpoint).l2Norm() /
                std::sqrt(permittivity_);
            if (neighbourFrequency > frequency &&
                neighbourFrequency < maxFrequency) {
                // Prevent double insertion by checking if the neighbour is
                // already in the queue. We explicitly check the
                // coordinates, as the frequency is a floating point number.
                bool found = false;
                for (auto iter = queue.begin(); iter != queue.end(); ++iter) {
                    if (std::get<1>(*iter) == n) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    queue.emplace(neighbourFrequency, n);
                }
            }
        }
    }
    // group the results
    return group(frequencies, 1e-5);
}

template <std::size_t DIM>
std::unique_ptr<typename BandStructure<DIM>::LineSet>
    HomogeneousBandStructure<DIM>::computeLines(
        LinearAlgebra::SmallVector<DIM> point1,
        LinearAlgebra::SmallVector<DIM> point2, double maxFrequency) const {
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

    for (LatticePoint<DIM> n :
         LatticePoint<DIM>::getNearestNeighbours(point1, reciprocalVectors_)) {
        considered.emplace(n);
        queue.emplace(n);
    }
    while (!queue.empty()) {
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
        if (fmin > maxFrequency) {
            continue;
        }
        // Add the point to the result
        LinearAlgebra::SmallVector<2> newMode({x, y});
        bool added = false;
        for (auto& mode : modes) {
            if ((mode.first - newMode).l2NormSquared() < 1e-10) {
                mode.second++;
                added = true;
            }
        }
        if (!added) {
            modes[newMode] = 1;
        }

        for (LatticePoint<DIM> n : p.getNeighbours()) {
            // Enqueue those neighbours that have not already been visited.
            if (considered.find(n) != considered.end()) {
                continue;
            }
            queue.push(n);
            considered.emplace(n);
        }
    }
    // construct the result
    std::unique_ptr<LineSet> result(new LineSet(l, permittivity_));
    for (auto& mode : modes) {
        result->addLine(mode.first[0], mode.first[1], mode.second);
    }
    return result;
}

template class HomogeneousBandStructure<2>;
template class HomogeneousBandStructure<3>;
