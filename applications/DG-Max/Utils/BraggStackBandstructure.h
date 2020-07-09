
#ifndef HPGEM_APP_BRAGGSTACKBANDSTRUCTURE_H
#define HPGEM_APP_BRAGGSTACKBANDSTRUCTURE_H

// So far only in 3D

#include "LinearAlgebra/SmallVector.h"
#include "BandStructure.h"

using namespace hpgem;


/// \brief Theoretical Bandstructure of a Bragg stack
///
/// Bandstructure of a Bragg stack using the theory from Yariv & Yeh (Optical
/// Waves in Crystals), section 6.2. There they use a stack in Z direction while
/// here we use a stack in X direction, hence some difference in variables.
/// Yariv and Yeh give a implicit condition for a valid TE mode in 6.2-24 of the
/// form cos(L k_x) - 0.5(A+D) = 0, where L is the period of the stack, k_x is
/// the wave vector in the X direction and A and D are two expressions from
/// 6.2-12, which depent on the transverse wavevector k_y and the angular
/// frequency omega. The equations for A and D change to those of 6.2-14 for the
/// TM modes.
///
/// For the implementation we note two properties:
///  - The bands are folded in the X-direction, as adding 2pi/L, the reciprocal
///    lattice vector, to k_z does not change cos(L k_x). We therefore do not
///    have to consider multiple lattice points in the X-direction, as all the
///    valid modes (omega, k_x, k_y) are already generated using the wavevector
///    k_x' that is shifted to the first Brillouin zone.
///  - The k_y vector is the transverse wavevector, which is assumed to be in
///    the Y-direction for the analysis. As the transverse direction can be
///    rotated freely this is completely valid. For implementation though we
///    need to fix the axis and need to support transverse wavevectors in the
///    whole transverse plane (k_y, k_z). This change is rather simple, by
///    replacing the k_y in the theory by sqrt(k_y^2 + k_z^2) (note that k_z
///    in Yariv & Yeh is what we here denote by k_x).
///  - Yariv & Yeh only use a periodic structure in the Bragg stack direction,
///    treating the transverse direction continuously. For computing band
///    structures we use a periodic structure in all three directions, we
///    therefore need to introduce a fictitious lattice in the transverse
///    direction. To keep this simple we use a cube as unit cell. This results
///    in some band folding.
///
class BraggStackBandstructure : public BandStructure<3> {
   public:
    BraggStackBandstructure(double eps1, double eps2, double fraction = 0.5);

    /// Compute the spectrum at a specific k-point in the first Brillouin zone.
    /// \param kpoint The point to compute the spectrum at
    /// \param omegaMax The maximum angular frequency to compute
    /// \return A map from frequency to to degeneracy of that frequency.
    std::map<double, std::size_t> computeSpectrum(
        LinearAlgebra::SmallVector<3> kpoint, double omegaMax) const final;

    /// \brief Set of Lines in the banstructure while traversing a line in
    /// k-space.
    class LineSet : public BandStructure<3>::LineSet {
       private:
        struct MultiMode;

       public:
        LineSet(const BraggStackBandstructure& structure, double l, double x1,
                double x2);

        double frequency(std::size_t line, double interpolation) const final;
        std::size_t multiplicity(std::size_t line) const final;
        std::size_t numberOfLines() const final;

        int lineType(std::size_t line) const final;
        std::string lineTitle(std::size_t line) const final;

        // TODO: Friend structure.
        void addLine(std::size_t multiplicity, double x, double y,
                     std::pair<std::size_t, std::size_t> modes);

       private:
        std::pair<const MultiMode*, std::size_t> findMode(
            std::size_t line) const;

        const BraggStackBandstructure& structure_;
        const double l_;
        const double x1_, x2_;
        std::vector<MultiMode> multiModes_;
        std::size_t totalNumberOfLines_;

        /// \brief A multi mode is a combination of modes sharing transverse
        /// wavevectors
        ///
        /// As discussed in the documentation of Bragg stack, the implicit
        /// definition of the valid modes, already contains the folding in the
        /// X direction. Hence, by fixing the wavevector (k_x, k_y, k_z) we
        /// implicitly define a series of modes (omega, k_x + 2npi/L).
        /// The transverse component of the line is defined by
        /// k_t = (k_y, k_z) = k - k_i, where k is the transverse component of
        /// the point in the band diagram and k_i is a reciprocal lattice point.
        /// As only the magnitude of the transverse component is part of the
        /// implicit definition of valid modes, there could be several
        /// reciprocal lattice points which have identical magnitude of the
        /// transverse wavevector when following the line of the LineSet (see
        /// for an example below). The resulting modes from these reciprocal
        /// lattice points would thus be degenerate.
        ///
        /// Thus a Multimode on a line through k-space describes:
        ///  - A set of modes TE and TM that have a fixed magnitude of the
        ///    transverse wavevector (that varies along the line) and have the
        ///    same reduced orthogonal wavevector k_x.
        ///  - The degeneracy or multiplicity of the lines due to multiple
        ///    reciprocal lattice vectors resulting in the same magnitude of the
        ///    transverse wavevector alogn the line
        ///
        /// To allow the handling of the degenerate modes we use a similar
        /// approach as with the HomogeneousBandStructure. When looking at the
        /// transverse wavevector from a point in reciprocal lattice space, we
        /// can write it in the form
        /// |k_t| = sqrt(y^2 + (x+tl)^2)
        /// Where l is the lenght of the line segment described by this lineset,
        /// t in [0,1] parameterizes where we are on the line and x,y are two
        /// values that describing the distance. Specifically, y is the distance
        /// between the reciprocal lattice point and its orthogonal projection
        /// on the (infinite) line, while x is the distance between the
        /// projected point and the first point on the line, with the sign such
        /// that x + l is the distance between the projected point and the
        /// second point on the line. Note that all these distances are in the
        /// transverse plane, disregarding any change in k_x along the line.
        ///
        /// As noted above, the equivalence between to modes hold when comparing
        /// the transverse wavevectors in magnitude. This is best explained by
        /// example, say we look at the LineSet between Gamma = [000] and
        /// X = [100]/2. The wavevector from the two reciprocal lattice points
        /// [010] and [0 -1 0] are [x -1 0] and [x 1 0] (where x between 0 and
        /// 0.5 as we move from Gamma to X). So the magnitude of the transverse
        /// component of the wavevector remains constant at 1 reciprocal
        /// lattice vector. If we look at the formula from Yariv and Yeh we only
        /// see the transverse wavevector in 6.2-6, where only its magnitude is
        /// used. Hence, the modes from [010] and [0 -1 0] along Gamma - X will
        /// overlap exactly. In fact, for a square lattice they will overlap
        /// with [001] and [0 0 -1] giving a four fold degeneracy. This will
        /// split when transitioning to the X - S=[110]/2 path. The transverse
        /// wavevector will decrease for reciprocal lattice point [010] and
        /// increase for [0 -1 0], points [001] and [0 0 -1] will keep the
        /// same distance. Hence we would the one MultiMode with fourfold
        /// degeneracy for Gamma-X splits in 3 MultiModes on X-S, two non
        /// degenerate modes and one with degeneracy 2.
        struct MultiMode {
            MultiMode(double x, double y, std::size_t teModes,
                      std::size_t tmModes, std::size_t multiplicity)
                : x_(x),
                  y_(y),
                  teModes_(teModes),
                  tmModes_(tmModes),
                  multiplicity_(multiplicity) {}

            double x_;
            double y_;
            std::size_t teModes_;
            std::size_t tmModes_;
            std::size_t multiplicity_;
        };
    };

    std::unique_ptr<typename BandStructure<3>::LineSet> computeLines(
        LinearAlgebra::SmallVector<3> point1,
        LinearAlgebra::SmallVector<3> point2, double maxFrequency) const final;

   private:
    const double eps1_;
    const double eps2_;
    /// \brief The fraction of the stack consisting of material with eps1.
    const double fraction_;

    // Function indicator function, which roots determine which paramaters
    // give a valid TE mode.
    double valuete(double omega, LinearAlgebra::SmallVector<3> k) const;
    // See valuete, but for TM modes.
    double valuetm(double omega, LinearAlgebra::SmallVector<3> k) const;

    /// \brief Find roots (angular frequencies) of the te or tm mode functions.
    ///
    /// Finds the roots of the functions defining valid TE/TM modes
    /// (valuete, valuetm), thereby findind the frequencies at which a TE/TM
    /// mode exists.
    ///
    /// \param k The k vector at which to find the roots
    /// \param omegamax The maximum angular frequency to consider
    /// \param tm Whether to look for roots of the TM (true) or TE (false)
    /// function.
    /// \param out The vector to which to add the found frequencies.
    void findRoots(LinearAlgebra::SmallVector<3> k, double omegamax, bool tm,
                   std::vector<double>& out) const;

    /// \brief Internal function for findRoots, where we look for the roots on
    ///  a small angular frequency interval.
    ///
    /// We assume that the interval is very small compared to the variation in
    /// valuete/tm and hence that there is at most one root or two very closely
    /// spaced roots.
    /// \param k The wavevector
    /// \param omin The minimum angular frequency
    /// \param omax The maximum angular frequency
    /// \param tm TM or TE mode?
    /// \param out Vector to output the results in.
    void findRootsInterval(LinearAlgebra::SmallVector<3> k, double omin,
                           double omax, bool tm,
                           std::vector<double>& out) const;

    /// \brief Similar to findRoots, but find the n-th root
    double findRoot(LinearAlgebra::SmallVector<3> k, std::size_t n,
                    bool tm) const;
};

#endif  // HPGEM_APP_BRAGGSTACKBANDSTRUCTURE_H
