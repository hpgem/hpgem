
#ifndef HPGEM_ANALYTICBANDSTRUCTURE_H
#define HPGEM_ANALYTICBANDSTRUCTURE_H

#include "LinearAlgebra/SmallVector.h"

#include <map>
#include <set>

/// \brief Computation of the analytic band structure of a homogeneous structure/
template<std::size_t DIM>
class HomogeneousBandStructure
{
public:
    ///
    /// \param reciprocalVectors The lattice vectors to use
    /// \param permittivity The Relative permitivity of the material
    HomogeneousBandStructure(std::array<LinearAlgebra::SmallVector<DIM>, DIM> reciprocalVectors, double permittivity = 1);

    /// \brief Compute the spectrum at a given point in k-space
    ///
    /// \param kpoint A point in the First Brillouin zone
    /// \param numberOfModes The number of modes to include (including multiplicity)
    /// \return A listing of the modes with their multiplicities
    std::map<double, std::size_t> computeSpectrum(LinearAlgebra::SmallVector<DIM> kpoint, int numberOfModes) const;

    /// \brief A dispersion line in a band structure between two points in k-space
    ///
    /// This is a line in the bandstructure diagram, with its multiplicity, as
    /// calculated on the line that passes two kpoints k1, k2. This restriction
    /// is needed to allow grouping the bands from multiple reciprocal lattice
    /// points in k-space that give rise to the higher multiplicity.
    class Line
    {
    public:
        Line(double l, double x, double y, double permittivity, std::size_t multiplicity)
            : l_(l), x_(x), y_(y), permittivity_(permittivity), multiplicity_ (multiplicity)
        {}

        Line(const Line& other) = default;
        Line(Line&& other) = default;

        /// Compute the frequency at the point interpolated between k1 and k2
        ///
        /// \param interpolation Interpolation constant (k = k1 + a(k2 - k1))
        /// \return The frequency of this band at the specified point.
        double frequency(double interpolation) const
        {
            double x0 = l_*interpolation + x_;
            return std::sqrt(y_*y_ +x0*x0) / std::sqrt(permittivity_);
        }
        /// \return Multiplicity without taking into account the polarization
        std::size_t multiplicity() const
        {
            return multiplicity_;
        }

        double getL() const
        {
            return l_;
        }

        double getX() const
        {
            return x_;
        }

        double getY() const
        {
            return y_;
        }

        double getPermittivity() const
        {
            return permittivity_;
        }


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
        double x_;
        double y_;
        double permittivity_;
        std::size_t multiplicity_;
    };

    std::vector<Line> computeLines(LinearAlgebra::SmallVector<DIM> point1, LinearAlgebra::SmallVector<DIM> point2, double maxFrequency) const;

private:
    const std::array<LinearAlgebra::SmallVector<DIM>, DIM> reciprocalVectors_;
    double permittivity_;
};


#endif //HPGEM_ANALYTICBANDSTRUCTURE_H
