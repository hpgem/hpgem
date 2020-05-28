
#ifndef UTILS_KSPACEPATH_h
#define UTILS_KSPACEPATH_h

#include "LinearAlgebra/SmallVector.h"

template <std::size_t DIM>
class KSpacePath {
   public:
    using KPoint = LinearAlgebra::SmallVector<DIM>;

    /// Construct a path in k-space from a set of points and the amount of steps
    /// used to walk between them.
    /// \param points The points to walk between
    /// \param steps The number of steps to take from each point to the next.
    KSpacePath(std::vector<KPoint> points, std::size_t steps);
    std::size_t totalNumberOfSteps() const;
    /// \brief The K vector at the indexed point.
    ///
    /// \param index The index of the point
    /// \return The associated k-vector
    KPoint k(std::size_t index) const;
    /// \brief The step between the indexed k vector and previous vector.
    ///
    /// This will be the value of k(index) - k(index-1), for index == 0 it will
    /// be k(0). \param index The index of the step. \return The step in k
    /// space.
    KPoint dk(std::size_t index) const;

    /// \brief Whether dk did change from index-1 to index.
    ///
    /// \param index The index for the change.
    /// \return Whether dk(index-1) != dk(index).
    bool dkDidChange(std::size_t index) const;

    /// \return The number of corner points in this k-space walk.
    std::size_t numberOfCornerPoints() const;

    /// \brief Corner point of this walk through k-space
    ///
    /// \param index The corner index.
    /// \return The indexed-corner point
    KPoint kcorner(std::size_t index) const;

    /// \brief Construct a path consisting of a single point.
    /// \param point The point to walk to.
    /// \return A KSpacePath with two points, the origin and the given point,
    /// reach in a single step.
    static KSpacePath<DIM> singleStepPath(KPoint point);
    /// \brief Create the standard high-symmetry point path for a cube.
    /// \param steps The number of steps in each direction.
    /// \param back Whether to return to the gamma point.
    /// \return The path for the cube
    static KSpacePath<DIM> cubePath(std::size_t steps, bool back = false);

   private:
    const std::vector<KPoint> points_;
    const std::size_t steps_;
};

#endif  // UTILS_KSPACEPATH_h
