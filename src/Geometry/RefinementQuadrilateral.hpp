#ifndef RefinementQuadrilateral_hpp
#define RefinementQuadrilateral_hpp

#include <string>
#include <vector>

#include "Geometry/ReferenceGeometry.hpp"
#include "Geometry/PhysicalGeometry.hpp"
#include "Geometry/PointPhysical.hpp"
#include "Geometry/RefinementGeometry.hpp"

namespace Geometry
{
    class RefinementQuadrilateral : public Geometry::RefinementGeometry<2>
    {
    public:
        typedef PointPhysical<2>                PointPhysicalT;
        typedef PhysicalGeometry<2>             PhysicalGeometryT;
        typedef ReferenceGeometry<2>            ReferenceGeometryT;
        typedef std::vector<PointPhysicalT>     VectorOfPointPhysicalsT;
        typedef std::vector<unsigned int>       VectorOfIndicesT;

        /// Constructors.
        RefinementQuadrilateral(const ReferenceGeometryT* const referenceGeometry,
                          const PhysicalGeometryT* const physicalGeometry)
                        : referenceGeometry_(referenceGeometry), physicalGeometry_(physicalGeometry)
        {
          std::cout << "RefinementQuadrilateral(referenceGeometry, physicalGeometry)\n";
        }

        RefinementQuadrilateral(const RefinementQuadrilateral& other)
                        : referenceGeometry_(other.referenceGeometry_), 
                          physicalGeometry_(other.physicalGeometry_)
        {
          std::cout << "RefinementQuadrilateral(other)\n";
        }

        virtual ~RefinementQuadrilateral()
        {}

        /// For debugging and checkpointing: a human-readable name.
        virtual std::string getName() const
        { 
          return "RefinementQuadrilateral";
        }

        //---------------------- Refinement definitions -----------------------------------------
        
        /// Number of new nodes due to a refinement.
        virtual unsigned int nrOfNewNodes(int refineType) const;
        
        /// Get all physical nodes: existing nodes and new nodes to be added.
        virtual void getAllNodes(int refineType, VectorOfPointPhysicalsT& nodes) const;
        
        /// Number of sub-elements due to the refinement.
        virtual unsigned int nrOfSubElements(int refineType) const;

        /// Assembly nodes for sub-element.
        virtual void subElementLocalNodeIndices(int refineType, unsigned int, VectorOfIndicesT& LocalNodeIdx) const;
        
        /// Local indices pairs of sub-elements connected by a sub-Internal Face.
        virtual void adjacentSubElementsPairs(int refineType,
                        VectorOfIndicesT& elemIdx1, VectorOfIndicesT& localFaceIdx1,
                        VectorOfIndicesT& elemIdx2, VectorOfIndicesT& localFaceIdx2) const;

        /// Number of sub-elements on a parent's face.
        virtual unsigned int nrOfSubElementsOnFace(int refineType, unsigned int faLocalIndex) const;

        /// Get sub-elements' local index on a parent's face.
        virtual void subElementsOnFace(int refineType, unsigned int faLocalIndex, VectorOfIndicesT& LocalNodeIdx) const;
        
        /// Get sub-face's local face number of on a parent's face.
        virtual unsigned int getLocalSubFaceNr(int refineType, unsigned int localFaceNr, unsigned int subElementIdx) const;

    private:
        RefinementQuadrilateral() : referenceGeometry_(NULL), physicalGeometry_(NULL) {}
        
        /// The physicalGeometry object contains pointers to the actual physical nodes, and
        /// a container of global node indexes.
        const PhysicalGeometryT*        physicalGeometry_;

        /// The corresponding referenceGeometry object
        const ReferenceGeometryT* const referenceGeometry_;
    };
}  // end of namespace Geometry

#endif