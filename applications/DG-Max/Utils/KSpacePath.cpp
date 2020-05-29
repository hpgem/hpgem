
#include "Logger.h"
#include "KSpacePath.h"

template <std::size_t DIM>
KSpacePath<DIM>::KSpacePath(std::vector<KSpacePath::KPoint> points,
                            std::size_t steps)
    : points_(points), steps_(steps) {}

template <std::size_t DIM>
std::size_t KSpacePath<DIM>::totalNumberOfSteps() const {
    // No steps are needed for the first point,
    return 1 + (points_.size() - 1) * steps_;
}

template <std::size_t DIM>
LinearAlgebra::SmallVector<DIM> KSpacePath<DIM>::k(std::size_t index) const {
    logger.assert_debug(index <= totalNumberOfSteps(),
                        "Index % is larger than the number of steps", index,
                        totalNumberOfSteps());
    if (index % steps_ == 0) {
        return points_[index / steps_];
    }
    // Linear interpolation
    std::size_t pointOffset = index / steps_;
    std::size_t stepOffset = index % steps_;
    double alpha = 1.0 * stepOffset / steps_;
    return alpha * points_[pointOffset + 1] +
           (1 - alpha) * points_[pointOffset];
}

template <std::size_t DIM>
LinearAlgebra::SmallVector<DIM> KSpacePath<DIM>::dk(std::size_t index) const {
    logger.assert_debug(index <= totalNumberOfSteps(),
                        "Index % is larger than the number of steps", index,
                        totalNumberOfSteps());
    if (index == 0) {
        logger(WARN, "Asking for dk at point 0");
        return points_[0];
    }
    std::size_t pointOffset = (index - 1) / steps_;
    return (points_[pointOffset + 1] - points_[pointOffset]) / steps_;
}

template <std::size_t DIM>
bool KSpacePath<DIM>::dkDidChange(std::size_t index) const {
    logger.assert_always(index < totalNumberOfSteps(),
                         "Index {} larger than total number of steps ", index,
                         totalNumberOfSteps());
    // When arriving at point a point for indices 0, steps_, 2steps_, etc.
    // so at 1, steps_ +1, 2steps_ + 1 we change direction. Additionally
    // for index == 0 dk is not valid, so we mark it as changed
    return index == 0 || (index - 1) % steps_ == 0;
}

template <std::size_t DIM>
std::size_t KSpacePath<DIM>::numberOfCornerPoints() const {
    return points_.size();
}

template <std::size_t DIM>
typename KSpacePath<DIM>::KPoint KSpacePath<DIM>::kcorner(
    std::size_t index) const {
    logger.assert_debug(index >= 0 && index < points_.size(), "Invalid index");
    return points_[index];
}

template <std::size_t DIM>
KSpacePath<DIM> KSpacePath<DIM>::singleStepPath(KSpacePath<DIM>::KPoint point) {
    std::vector<KSpacePath<DIM>::KPoint> points(1);
    points[0] = point;
    KSpacePath<DIM> path(points, 1);
    return path;
}

// TODO: This does not really seem nice.
template <std::size_t DIM>
KSpacePath<DIM> KSpacePath<DIM>::cubePath(std::size_t steps, bool back) {
    std::vector<KSpacePath<DIM>::KPoint> points(DIM + (back ? 2 : 1));
    KSpacePath<DIM>::KPoint point;
    // 0
    points[0] = point;
    if (back) {
        points[DIM + 1] = point;
    }
    // (pi, 0, 0)
    point[0] = M_PI;
    points[1] = point;
    // (pi, pi, 0)
    point[1] = M_PI;
    points[2] = point;
    if (DIM == 3) {
        // (pi, pi, pi)
        point[2] = M_PI;
        points[3] = point;
    }
    // TODO: Apparently this is incorrect
    KSpacePath<DIM> path(points, steps);
    return path;
}

// Explicitly instantiate for
template class KSpacePath<2>;
template class KSpacePath<3>;
