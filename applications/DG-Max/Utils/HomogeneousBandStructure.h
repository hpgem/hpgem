
#ifndef HPGEM_ANALYTICBANDSTRUCTURE_H
#define HPGEM_ANALYTICBANDSTRUCTURE_H

#include "LinearAlgebra/SmallVector.h"

#include <map>

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
    std::map<double, std::size_t> computeSpectrum(LinearAlgebra::SmallVector<DIM> kpoint, int numberOfModes);

private:
    std::array<LinearAlgebra::SmallVector<DIM>, DIM> reciprocalVectors_;
    double permittivity_;
};


#endif //HPGEM_ANALYTICBANDSTRUCTURE_H
