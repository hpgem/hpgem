#include "DivDGMaxEVTestCase.h"

#include "DGMaxLogger.h"
#include "DGMaxProgramUtils.h"

namespace DGMax {
template <std::size_t DIM>
std::unique_ptr<AbstractEigenvalueResult<DIM> >
    DivDGMaxEVTestCase<DIM>::runInternal(std::size_t level) {

    logger.assert_always(level < meshFileNames_.size(), "No such mesh");

    Base::ConfigurationData configData(2, 1);
    auto mesh = DGMax::readMesh<DIM>(
        meshFileNames_[level], &configData,
        [&](const Geometry::PointPhysical<DIM> &p) {
            return jelmerStructure(p, testCase_.getStructureId());
        });
    DGMaxLogger(INFO, "Loaded mesh % with % local elements.",
                meshFileNames_[level], mesh->getNumberOfElements());
    KSpacePath<DIM> path =
        KSpacePath<DIM>::singleStepPath(testCase_.getKPoint());
    EigenValueProblem<DIM> input(path, testCase_.getNumberOfEigenvalues());

    DivDGMaxEigenvalue<DIM> solver(*mesh, this->order_, this->stab_);
    return solver.solve(input);
}

// Template instantiation
template class DivDGMaxEVTestCase<2>;
template class DivDGMaxEVTestCase<3>;

}  // namespace DGMax