
#ifndef HPGEM_ANALYTICBANDSTRUCTURE_H
#define HPGEM_ANALYTICBANDSTRUCTURE_H

#include "LinearAlgebra/SmallVector.h"
#include "BandStructure.h"

#include <set>

/// \brief Computation of the analytic band structure of a homogeneous
/// structure/
template <std::size_t DIM>
class HomogeneousBandStructure : public BandStructure<DIM> {
   public:
    ///
    /// \param reciprocalVectors The lattice vectors to use
    /// \param permittivity The Relative permitivity of the material
    HomogeneousBandStructure(
        std::array<LinearAlgebra::SmallVector<DIM>, DIM> reciprocalVectors,
        double permittivity = 1);

    /// \brief Compute the spectrum at a given point in k space
    ///
    /// \param kpoint The k-point
    /// \param omegaMax  The maximum frequency
    /// \return The frequencies in increasing order at the k-point up to
    /// omegaMax.
    std::vector<double> computeLinearSpectrum(
        LinearAlgebra::SmallVector<DIM> kpoint, double omegaMax) const final;

    /// \brief Compute the grouped spectrum at a given point in k-space
    ///
    /// \param kpoint A point in the First Brillouin zone
    /// \param maxFrequency The maximum frequency of the mode
    /// \return A listing of the modes with their multiplicities
    std::map<double, std::size_t> computeSpectrum(
        LinearAlgebra::SmallVector<DIM> kpoint,
        double maxFrequency) const final;

    /// \brief A set of dispersion lines in a band structure between two points
    /// in k-space
    ///
    /// This is the set of lines in the bandstructure (including degeneracy) on
    /// a line that passes through two points in kspace k1, k2. This latter
    /// restriction allows for grouping of the bands from multiple reciprocal
    /// lattice points together, when they give the same band.
    class LineSet : public BandStructure<DIM>::LineSet {
       public:
        LineSet(double l, double permittivity)
            : l_(l), permittivity_(permittivity) {}

        LineSet(const LineSet& other) = default;
        LineSet(LineSet&& other) = default;

        /// Compute the frequency at the point interpolated between k1 and k2
        ///
        /// \param line The line number
        /// \param interpolation Interpolation constant (k = k1 + a(k2 - k1))
        /// \return The frequency of this band at the specified point.
        double frequency(std::size_t line, double interpolation) const final {
            double x0 = l_ * interpolation + xs_[line];
            double y = ys_[line];
            return std::sqrt(y * y + x0 * x0) / std::sqrt(permittivity_);
        }
        /// \return Multiplicity without taking into account the polarization
        std::size_t multiplicity(std::size_t line) const final {
            return multiplicities_[line];
        }

        double getPermittivity() const { return permittivity_; }

        void addLine(double x, double y, std::size_t multiplicity) {
            xs_.emplace_back(x);
            ys_.emplace_back(y);
            multiplicities_.emplace_back(multiplicity);
        }

        std::size_t numberOfLines() const final { return xs_.size(); }

       private:
        // The frequency omega for a band at point k is defined as
        // omega = |k - k_l| c/n
        // where c the speed of light (assumed to be 1), n the refractive index
        // of the material (sqrt(relative permittivity) and k_l is the lattice
        // point from which this band originates.
        //
        // When we restrict this to a line in k-space between two points k1 and
        // k2, we can expand the length by projecting the point k_l on the line
        // through k1 and k2. Given this projection point kp the distance is
        // given by |k - k_l|^2 = |k - kp|^2 + |kp - k_l|^2 (pythagoras).
        // We can directly see that |kp - kl| is constant independent on where
        // the point k is situated on the line, while |kp - kl| varies.
        // Moreover, there might be multiple lattice points k_l with the same
        // value for |kp - k_l|, resulting in a band with higher multiplicity.
        //
        //
        // If we look at the plane intersecting the three points k1,k2 and k_l,
        // we can represent the situation as follows
        //
        //  kp ---- k1 --[k somewhere]--k2
        //   |
        //   |
        //   |
        //   |
        //  k_l
        //
        // Point kp can easily be found, let dk = k1 - k_l and
        // kdir = k2 - k1 (normalized)
        // then
        // kp = k_l + dk - <dk, kdir>kdir = k1 - <dk, kdir>kdir
        //
        // Hence the second norm |kp - k_l| = |dk - <dk, kdir>kdir|
        // For the first norm we use the linear interpolation of
        // k = k1 + a(k2 - k1) and thus
        // k - kp = a(k2 - k1) + k1 - kp.
        // Taking the inner product with the unit vector kdir (collinear with
        // k2 - k1 and k1 - kp) we get
        // |k - kp| = |a<k2 - k1, kdir> + <k1 - kp, kdir>|
        // Note now that <k2 - k1, kdir> is the length of the line between
        // k1 and k2 and thus fixed, and
        // <k1 - kp, kdir> = <<dk, kdir>kdir, kdir> = <dk, kdir>
        // Where we used that kdir is a unit vector.
        // Hence we have
        // |k - k_l|^2 = |dk - <dk, kdir>kdir|^2 + |al + <dk, kdir>|^2
        // with l = |k2 - k1|
        // So instead of saving the point k_l, which would make this line object
        // dependent of the lattice point, and thus prevent having higher
        // multiplicity, we instead store three numbers
        // l = |k2 - k1|
        // x = <dk, kdir>
        // y = |dk - <dk, kdir>kdir|
        // and the previous diagram can thus be expanded to
        //   ---x---> --------l--------->
        //  kp ---- k1 --[k somewhere]--k2
        // ^ |    # -kdir->
        // | |   /
        // y |  dk
        // | | /
        //  k_l
        double l_;
        double permittivity_;
        // TODO: This should be struct
        std::vector<double> xs_;
        std::vector<double> ys_;
        std::vector<std::size_t> multiplicities_;
    };

    std::unique_ptr<typename BandStructure<DIM>::LineSet> computeLines(
        LinearAlgebra::SmallVector<DIM> point1,
        LinearAlgebra::SmallVector<DIM> point2,
        double maxFrequency) const final;

   private:
    const std::array<LinearAlgebra::SmallVector<DIM>, DIM> reciprocalVectors_;
    double permittivity_;
};

#endif  // HPGEM_ANALYTICBANDSTRUCTURE_H
