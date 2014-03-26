#ifndef RefinementPyramid_hpp
#define RefinementPyramid_hpp

#include "Geometry/ReferenceGeometry.hpp"
#include "Geometry/PhysicalGeometry.hpp"
#include "Geometry/PointReference.hpp"
#include "Geometry/PointPhysical.hpp"
#include "LinearAlgebra/Matrix.hpp"
#include "Geometry/RefinementGeometry.hpp"

namespace Geometry
{
    class RefinementPyramid : public Geometry::RefinementGeometry
    {
    public:
        typedef unsigned int                    DimT;
        typedef Point                        PointT;
        typedef PointPhysical                PointPhysicalT;
        typedef PointReference               PointReferenceT;
        typedef PhysicalGeometry             PhysicalGeometryT;
        typedef ReferenceGeometry            ReferenceGeometryT;
        typedef std::vector<PointPhysicalT>     VectorOfPointPhysicalsT;
        typedef std::vector<unsigned int>       VectorOfIndicesT;

        /// Constructors.
        RefinementPyramid(const ReferenceGeometryT* const referenceGeometry,
                          const PhysicalGeometryT* const physicalGeometry)
                        : referenceGeometry_(referenceGeometry), physicalGeometry_(physicalGeometry)
        {
          std::cout << "RefinementPyramid(referenceGeometry, physicalGeometry)\n";
        }

        RefinementPyramid(const RefinementPyramid& other)
                        : referenceGeometry_(other.referenceGeometry_), 
                          physicalGeometry_(other.physicalGeometry_)
        {
          std::cout << "RefinementPyramid(other)\n";
        }

        virtual ~RefinementPyramid()
        {}

        /// For debugging and checkpointing: a human-readable name.
        virtual std::string getName() const
        { 
          return "RefinementPyramid";
        }

        //---------------------- Refinement definitions -----------------------------------------
        
        /// Number of new nodes due to the refinement that should be added to
        /// the vector of localNodeIndices
        virtual DimT nrOfNewNodes(int refineType) const { return 0; }
        
        /// New physical nodes due to refinement to the nodes vector
        /// \param newPoints On input, this vector contains the element's physical nodes.  
        /// On exit, the new physical nodes are added.
        virtual void getAllNodes(int refineType, VectorOfPointPhysicalsT& nodes) const {}
        
        /// Number of sub-elements due to the refinement
        virtual DimT nrOfSubElements(int refineType) const { return 0; }

        /// Assembly nodes for sub-element
        virtual void subElementLocalNodeIndices(int refineType, DimT iSubElement, VectorOfIndicesT& LocalNodeIdx) const {}
        
        /// Local indices pairs of sub-elements connected by a sub-Internal Face
        virtual void adjacentSubElementsPairs(int refineType,
                        VectorOfIndicesT& elemIdx1, VectorOfIndicesT& localFaceIdx1,
                        VectorOfIndicesT& elemIdx2, VectorOfIndicesT& localFaceIdx2) const {}

        /// Number of sub-elements on a parent's face.
        virtual DimT nrOfSubElementsOnFace(int refineType, DimT faLocalIndex) const { return 0; }

        /// Get sub-elements' local index on a parent's face.
        virtual void subElementsOnFace(int refineType, DimT faLocalIndex, VectorOfIndicesT& localSubElemIdx) const {}
        
        /// Get sub-face's local face number of on a parent's face.
        virtual DimT getLocalSubFaceNr(int refineType, DimT localFaceNr, DimT subElementIdx) const { return 0; }

    private:
        RefinementPyramid() : referenceGeometry_(NULL), physicalGeometry_(NULL) {}

        /// The physicalGeometry object contains pointers to the actual physical nodes, and
        /// a container of global node indexes.
        const PhysicalGeometryT*        physicalGeometry_;

        /// The corresponding referenceGeometry object
        const ReferenceGeometryT* const referenceGeometry_;
    };
}
#endif
