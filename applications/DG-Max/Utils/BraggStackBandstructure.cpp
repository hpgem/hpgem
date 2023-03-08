
#include "BraggStackBandstructure.h"

#include "BandstructureUtils.h"

#include <set>
#include <queue>

using namespace hpgem;

template <std::size_t DIM>
BraggStackBandstructure<DIM>::BraggStackBandstructure(double eps1, double eps2,
                                                      double fraction,
                                                      double transverseSize)
    : eps1_(eps1),
      eps2_(eps2),
      fraction_(fraction),
      computeTE_(true),
      computeTM_(true),
      transverseSize_(transverseSize) {}

template <std::size_t DIM>
std::vector<double> BraggStackBandstructure<DIM>::computeLinearSpectrum(
    LinearAlgebra::SmallVector<DIM> kpoint, double omegaMax) const {
    for (std::size_t i = 0; i < DIM; ++i)
        logger.assert_always(std::abs(kpoint[i]) <= (M_PI + 1e-8),
                             "Only implemented for k in first Brillouin zone.");
    // Start with only for this specific kpoint.
    std::vector<double> freqs;
    std::set<LatticePoint<DIM - 1>> considered;
    std::queue<LatticePoint<DIM - 1>> toProcess;

    LatticePoint<DIM - 1> origin{};
    considered.emplace(origin);
    toProcess.emplace(origin);

    double kp = kpoint[0];
    LinearAlgebra::SmallVector<DIM - 1> ktransverse;
    for (std::size_t i = 1; i < DIM; ++i) {
        ktransverse[i - 1] = kpoint[i];
    }

    while (!toProcess.empty()) {
        LatticePoint<DIM - 1> offset = toProcess.front();
        toProcess.pop();
        LinearAlgebra::SmallVector<DIM - 1> ktransverseRel = ktransverse;
        for (std::size_t i = 0; i < DIM - 1; ++i) {
            ktransverseRel[i] += 2 * M_PI * offset[i] / transverseSize_;
        }
        double ktransverseRelMag = ktransverseRel.l2Norm();
        std::size_t size = freqs.size();
        if (computeTE_) {
            findRoots(kp, ktransverseRelMag, omegaMax, false, freqs);
        }
        if (computeTM_) {
            findRoots(kp, ktransverseRelMag, omegaMax, true, freqs);
        }
        if (freqs.size() != size) {
            // This offset contributed, add its non considered neighbours.
            for (auto neighbour : offset.getNeighbours()) {
                if (considered.find(neighbour) == considered.end()) {
                    considered.insert(neighbour);
                    toProcess.push(neighbour);
                }
            }
        }
    }
    std::sort(freqs.begin(), freqs.end());
    return freqs;
}

template <std::size_t DIM>
std::map<double, std::size_t> BraggStackBandstructure<DIM>::computeSpectrum(
    LinearAlgebra::SmallVector<DIM> kpoint, double omegaMax) const {
    std::vector<double> frequencies = computeLinearSpectrum(kpoint, omegaMax);
    return group(frequencies, 1e-12);
}

/// \brief Simple bisection algorithm to find the root of a functions
///
/// Bisection to find the root of a function, it is expected that the function
/// differs in sign between the two endpoints of the starting interval
///
/// \param function The function to find the root of
/// \param left Left end of the interval
/// \param right Right end of the interval
/// \param tol Tolerance in the coordinates
/// \return Coordinate of the root
double bisect(std::function<double(double)> function, double left, double right,
              double tol = 1e-12) {
    double fleft = function(left);
    double fright = function(right);
    while (std::abs(right - left) > tol) {
        double mid = (left + right) / 2;
        double fmid = function(mid);
        if (fmid * fleft < 0) {
            right = mid;
            fright = fmid;
        } else if (fmid * fright < 0) {
            left = mid;
            fleft = fmid;
        } else {
            return mid;
        }
    }
    return (left + right) / 2;
}

template <std::size_t DIM>
void BraggStackBandstructure<DIM>::findRoots(double kp, double kt,
                                             double omegamax, bool tm,
                                             std::vector<double>& out) const {
    const std::size_t POINTS = 201;
    const double domega = omegamax / (POINTS - 1);
    for (std::size_t i = 0; i < POINTS - 1; ++i) {
        double o1 = i * domega;
        double o2 = o1 + domega;
        findRootsInterval(kp, kt, o1, o2, tm, out);
    }
    // Return via parameter
}

template <std::size_t DIM>
double BraggStackBandstructure<DIM>::findRoot(double kp, double kt,
                                              std::size_t n, bool tm) const {
    const double STEP_SIZE = 0.025;
    std::vector<double> roots;
    double omin = 0;
    while (roots.size() <= n) {
        double omax = omin + STEP_SIZE;
        findRootsInterval(kp, kt, omin, omax, tm, roots);
        omin = omax;
    }
    return roots[n];
}

template <std::size_t DIM>
void BraggStackBandstructure<DIM>::findRootsInterval(
    double kp, double kt, double omin, double omax, bool tm,
    std::vector<double>& out) const {
    auto value = [=](double o) {
        return tm ? valuetm(o, kp, kt) : valuete(o, kp, kt);
    };
    double fval1 = value(omin);
    double fval2 = value(omax);
    if (std::abs(fval1) < 1e-12) {
        // Left endpoint is small enough to expect there to be a root very
        // nearby either due to finite precision in computing f, or in the point
        // omin. We do not do the same check for fval2 to prevent duplicate
        // detection.
        out.emplace_back(omin);
    } else if (std::abs(fval2) > 1e-8 && fval1 * fval2 < 0) {
        // We have different signs of fval1 and fval2 thus there must be a zero
        // crossing inbetween omin and omax.
        out.emplace_back(bisect(value, omin, omax));
    } else if (std::abs(fval2) < 0.5 && std::abs(fval1) < 0.5) {
        // Two small values of f, it might be only touching or barely crossing
        // zero. We use golden section search to find the minimum of f on the
        // interval. If in the process of doing so we detect a zero crossing,
        // we switch back to bisection to find those locations.

        // multiply the values of valuete/valuetm by the original sigh of fval,
        // this ensures that we are always having positive values and are
        // looking for a minimum.
        double sign = fval1 > 0 ? 1 : -1;
        std::function<double(double)> func = [=](double o) {
            return sign * value(o);
        };

        // End points of the interval
        double a = omin, b = omax;
        double fa = sign * fval1, fb = sign * fval2;

        // From the fact that fval1 and fval2 are small, we do not know if there
        // is a minimum on the interval, or maybe an minimum just outside the
        // interval. Therefore, we compare the values just outside of the
        // interval, for a minimum to be in the interval these should be larger
        // than the endpoints. This is similar to computing the derivatives at
        // the endpoints and observe that it is negative for the left endpoint
        // of the interval and positive for the right one.
        double tmin = omin - 0.25 * (omax - omin);
        double tmax = omax + 0.25 * (omax - omin);
        double fmin = func(tmin);
        double fmax = func(tmax);
        if (func(tmin) < fa || func(tmax) < fb) {
            return;
        }

        // Golden ratio
        const double phi = (1 + std::sqrt(5)) / 2;
        // Two internal points that are used to determine how to shrink the
        // interval
        double c = b - (b - a) / phi;
        double d = a + (b - a) / phi;
        double fc = func(c), fd = func(d);

        bool rootsFound = false;
        // Golden section loop
        while (std::abs(c - d) > 1e-12 && !rootsFound) {
            if (fc < 0) {
                // Negative value of f found at c while f(a) and f(b) are
                // positive. Hence, we must have two positions d,e with
                // a < d < c < e < b such that f(d) = f(e) = 0.
                out.emplace_back(bisect(func, a, c));
                out.emplace_back(bisect(func, c, b));
                rootsFound = true;
            } else  if (fd < 0) {
                // Similar to the case fc < 0, but now with fd.
                out.emplace_back(bisect(func, a, d));
                out.emplace_back(bisect(func, d, b));
                rootsFound = true;
            } else if (fc < fd) {
                // Actual golden section case 1:
                // fc < fd means that the minimum must be between a < c < d
                // and point b (>d) can be eliminated as candidate
                b = d;
                fb = fd;
                d = c;
                fd = fc;
                c = b - (b - a) / phi;
                fc = func(c);
            } else {
                // Golden section case 2: Point a can be discarded
                a = c;
                fa = fc;
                c = d;
                fc = fd;
                d = a + (b - a) / phi;
                fd = func(d);

            }
        }
        if (!rootsFound && std::abs(fa) < 1e-12) {
            // fa is small enough to assume that there is an x very near to a
            // such that f(x) = 0.
            double mid = (a + b) / 2;
            out.emplace_back(mid);
            out.emplace_back(mid);
        }
    }
}

template <std::size_t DIM>
std::unique_ptr<typename BandStructure<DIM>::LineSet>
    BraggStackBandstructure<DIM>::computeLines(
        LinearAlgebra::SmallVector<DIM> point1,
        LinearAlgebra::SmallVector<DIM> point2, double maxFrequency) const {
    // Compute the lines in similar way to the Homogeneous case, but here we
    // only have a fictitious lattice in 2 directions. Moreover we have to
    // consider that on a transverse wavevector a whole set of valid modes is
    // constructed by varying kx (the orthogonal wavevector).

    for (std::size_t i = 0; i < DIM; ++i) {
        logger.assert_always(std::abs(point1[i]) <= (M_PI + 1e-12),
                             "Only implemented for k in first Brillouin zone.");
        logger.assert_always(std::abs(point2[i]) <= (M_PI + 1e-12),
                             "Only implemented for k in first Brillouin zone.");
    }

    LinearAlgebra::SmallVector<DIM - 1> tpoint1, tpoint2, kdir;
    // Extract the transverse components of point1, point2
    for (std::size_t i = 1; i < DIM; ++i) {
        tpoint1[i - 1] = point1[i];
        tpoint2[i - 1] = point2[i];
    }
    // Note that unlike the vacumm case, the line Gamma-X has no change in the
    // transverse direction. Hence, l2 and kdir are both zero.
    kdir = tpoint2 - tpoint1;
    const double l2 = kdir.l2Norm();
    if (std::abs(l2) > 1e-12) {
        kdir /= l2;
    }

    // Set of points that are considered (either in the queue or were in the
    // queue at some previous point).
    std::set<LatticePoint<DIM - 1>> considered;
    // Queue of points that still have to be processed.
    std::queue<LatticePoint<DIM - 1>> queue;
    // Identify the resulting lines/modes by the two coordinates, x,y
    // and associate with it the multiplicity.
    std::map<LinearAlgebra::SmallVector<2>, std::size_t> modes;
    std::map<LinearAlgebra::SmallVector<2>, std::pair<std::size_t, std::size_t>>
        modes2;

    LatticePoint<DIM - 1> start;
    considered.emplace(start);
    queue.emplace(start);

    while (!queue.empty()) {
        LatticePoint<DIM - 1> p = queue.front();
        queue.pop();
        // Reciprocal wavevector in the transverse direction.
        LinearAlgebra::SmallVector<DIM - 1> dk1 = tpoint1;
        LinearAlgebra::SmallVector<DIM - 1> dk2 = tpoint2;
        for (std::size_t i = 0; i < DIM - 1; ++i) {
            dk1[i] -= 2 * M_PI * p[i] / transverseSize_;
            dk2[i] -= 2 * M_PI * p[i] / transverseSize_;
        }
        // Distance along the line
        double x = dk1 * kdir;
        // Distance to the line
        double y = (dk1 - x * kdir).l2Norm();
        double xmin = intervalDist(x, -l2, 0);
        double mindist = std::sqrt(xmin * xmin + y * y);

        // Unlike the with the homogeneous bandstructure we can not directly
        // compute the frequency of the contribution of a wavevector because:
        //  - There are many (with different kz)
        //  - The frequency is only defined implicitly.
        // Hence we use a rather crude estimate here to exclude those
        // wavevectors which will definitely not have a valid frequency below
        // the maximum frequency.
        double k1z = -(eps1_ * maxFrequency * maxFrequency - mindist * mindist);
        double k2z = -(eps2_ * maxFrequency * maxFrequency - mindist * mindist);
        // Treshold for a TE mode to exist, as tresh <= 0.5(A+D) = cos(kz);
        double tresh = std::cosh(std::sqrt(k1z) * fraction_ +
                                 std::sqrt(k2z) * (1 - fraction_));
        if (k1z >= 0 && k2z >= 0 && std::abs(tresh) > 1) {
            continue;
        }

        // Add the wavevector to the result
        LinearAlgebra::SmallVector<2> newMode({x, y});
        bool added = false;
        for (auto& mode : modes) {
            if ((mode.first - newMode).l2NormSquared() < 1e-12) {
                mode.second++;
                added = true;
            }
        }
        if (!added) {
            // It is a new wavevector/reciprocal lattice point. We have to
            // determine how many modes have to be included in the banddiagram,
            // or at least an upper bound to this. For this we take the maximum
            // of the amount of bands at either end of the line.
            modes[newMode] = 1;  // Single degeneracy
            std::vector<double> roots;
            // Transverse wavevector for k-point1 relative to lattice point kp
            double kt_rel1 = dk1.l2Norm();
            if (computeTE_) {
                findRoots(point1[0], kt_rel1, maxFrequency, false, roots);
            }
            std::size_t teRoots = roots.size();
            roots.clear();
            if (computeTM_) {
                findRoots(point1[0], kt_rel1, maxFrequency, true, roots);
            }
            std::size_t tmRoots = roots.size();
            roots.clear();
            // Transverse wavevector for k-point2 relative to lattice point kp
            double kt_rel2 = dk2.l2Norm();
            if (computeTE_) {
                findRoots(point2[0], kt_rel2, maxFrequency, false, roots);
            }
            teRoots = std::max(teRoots, roots.size());
            roots.clear();
            if (computeTM_) {
                findRoots(point2[0], kt_rel2, maxFrequency, true, roots);
            }
            tmRoots = std::max(tmRoots, roots.size());
            // If ky = kz = 0 then the TE and TM modes will overlap.
            // We test if ky and kz of the line don't change by checking l2
            // We test if the lattice point is on the line by checking y
            if (std::abs(l2) > 1e-12 || std::abs(y) > 1e-12 ||
                !(computeTM_ && computeTE_)) {
                modes2[newMode] = std::make_pair(teRoots, tmRoots);
            } else {
                // Thus remove tmMode but add 1 to multiplicity
                if (computeTE_ && computeTM_) {
                    modes2[newMode] = std::make_pair(teRoots, 0);
                    modes[newMode]++;
                }
            }
        }

        // Enqueue those neighbours that have not already been visited.
        for (LatticePoint<DIM - 1> n : p.getNeighbours()) {
            if (considered.find(n) != considered.end()) {
                continue;
            }
            queue.push(n);
            considered.emplace(n);
        }
    }

    // Compute line set.
    std::unique_ptr<BraggStackBandstructure<DIM>::LineSet> result(
        new LineSet(*this, l2, point1[0], point2[0]));
    for (auto& mode : modes) {
        std::pair<std::size_t, std::size_t> ms = modes2[mode.first];
        result->addLine(mode.second, mode.first[0], mode.first[1], ms);
    }
    return result;
}

template <std::size_t DIM>
double BraggStackBandstructure<DIM>::valuetm(double omega, double kp,
                                             double kt) const {
    // See Yarev & Yeh 6.2-24 for the original formula
    // + 6.2-12 for the coefficients
    // Note they use the opposite TE/TM convention from photonics
    double kt2 = kt * kt;
    std::complex<double> k1z =
        std::sqrt(std::complex<double>(eps1_ * omega * omega - kt2));
    std::complex<double> k2z =
        std::sqrt(std::complex<double>(eps2_ * omega * omega - kt2));
    std::complex<double> kk = k1z / k2z + k2z / k1z;

    if (omega < 1e-6 && kt2 < 1e-6)  // TODO: Better bounds
    {
        kk = std::sqrt(eps1_ / eps2_) + std::sqrt(eps2_ / eps1_);
    }

    double b = 1 - fraction_;
    // (A + D)/2
    double ad2 = (std::cos(k1z * fraction_) * std::cos(k2z * b) -
                  0.5 * kk * std::sin(k1z * fraction_) * std::sin(k2z * b))
                     .real();
    return std::cos(kp) - ad2;
}

template <std::size_t DIM>
double BraggStackBandstructure<DIM>::valuete(double omega, double kp,
                                             double kt) const {
    // See Yarev & Yeh 6.2-24 for the original formula
    // + 6.2-14 for the coefficients
    // Note they use the opposite TE/TM convention from photonics
    double kt2 = kt * kt;
    std::complex<double> k1z =
        std::sqrt(std::complex<double>(eps1_ * omega * omega - kt2));
    std::complex<double> k2z =
        std::sqrt(std::complex<double>(eps2_ * omega * omega - kt2));
    std::complex<double> kk =
        (eps2_ * k1z) / (eps1_ * k2z) + (eps1_ * k2z) / (eps2_ * k1z);
    if (omega < 1e-6 && kt2 < 1e-6)  // TODO: Better bounds
    {
        kk = std::sqrt(eps1_ / eps2_) + std::sqrt(eps2_ / eps1_);
    }

    double b = 1 - fraction_;
    double ad2 = (std::cos(k1z * fraction_) * std::cos(k2z * b) -
                  0.5 * kk * std::sin(k1z * fraction_) * std::sin(k2z * b))
                     .real();
    // (A + D)/2
    return std::cos(kp) - ad2;
}

// Lineset
template <std::size_t DIM>
BraggStackBandstructure<DIM>::LineSet::LineSet(
    const BraggStackBandstructure<DIM>& structure, double l, double x1,
    double x2)
    : structure_(structure), l_(l), x1_(x1), x2_(x2), totalNumberOfLines_(0) {}

template <std::size_t DIM>
double BraggStackBandstructure<DIM>::LineSet::frequency(
    std::size_t line, double interpolation) const {
    logger.assert_debug(line >= 0 && line <= numberOfLines(),
                        "Invalid line number");
    // Find line number
    auto selection = findMode(line);
    const MultiMode& multiMode = *(selection.first);
    std::size_t index = selection.second;

    std::size_t teModes = multiMode.teModes_;
    // Found the major mode
    bool tm = index >= teModes;
    if (tm) {
        index -= teModes;
    }
    std::vector<double> roots;
    // Similar to the Homogeneous case, compute the wavevector distance in the
    // transverse plane.
    double x = multiMode.x_, y = multiMode.y_;
    x += interpolation * l_;
    double z = x1_ + interpolation * (x2_ - x1_);
    LinearAlgebra::SmallVector<3> k({z, x, y});
    return structure_.findRoot(z, std::sqrt(x * x + y * y), index, tm);
}

template <std::size_t DIM>
void BraggStackBandstructure<DIM>::LineSet::addLine(
    std::size_t multiplicity, double x, double y,
    std::pair<std::size_t, std::size_t> modes) {
    multiModes_.emplace_back(x, y, modes.first, modes.second, multiplicity);
    totalNumberOfLines_ += modes.first + modes.second;
}

template <std::size_t DIM>
std::size_t BraggStackBandstructure<DIM>::LineSet::multiplicity(
    std::size_t line) const {
    return findMode(line).first->multiplicity_;
}

template <std::size_t DIM>
std::size_t BraggStackBandstructure<DIM>::LineSet::numberOfLines() const {
    return totalNumberOfLines_;
}

template <std::size_t DIM>
int BraggStackBandstructure<DIM>::LineSet::lineType(std::size_t line) const {
    std::pair<const MultiMode*, std::size_t> mode = findMode(line);
    if (std::abs(mode.first->x_) < 1e-8 && std::abs(mode.first->y_) < 1e-8 &&
        std::abs(l_) < 1e-8) {
        // TE+TM modes overlap when the transverse component of the wavevector
        // is zero.
        return 3;
    }
    if (mode.second < mode.first->teModes_) {
        // TE Mode
        return 1;
    }
    // TM Mode
    return 2;
}

template <std::size_t DIM>
std::string BraggStackBandstructure<DIM>::LineSet::lineTitle(
    std::size_t line) const {
    std::ostringstream out;
    int type = lineType(line);
    switch (type) {
        case 1:
            out << "TE ";
            break;
        case 2:
            out << "TM ";
            break;
        case 3:
            out << "TE+TM ";
            break;
        default:
            out << "Unknown";
            break;
    }
    std::pair<const MultiMode*, std::size_t> mode = findMode(line);
    out << mode.first->multiplicity_;
    return out.str();
}

template <std::size_t DIM>
std::pair<const typename BraggStackBandstructure<DIM>::LineSet::MultiMode*,
          std::size_t>
    BraggStackBandstructure<DIM>::LineSet::findMode(std::size_t line) const {
    logger.assert_debug(line >= 0 && line <= numberOfLines(),
                        "Invalid line number");
    for (const auto& multiMode : multiModes_) {
        if (line < multiMode.teModes_ + multiMode.tmModes_) {
            return std::make_pair(&multiMode, line);
        }
        line -= multiMode.teModes_ + multiMode.tmModes_;
    }
    logger.assert_always(false, "Invalid internal structure");
    return std::make_pair(nullptr, -1);
}

template class BraggStackBandstructure<2>;
template class BraggStackBandstructure<3>;
