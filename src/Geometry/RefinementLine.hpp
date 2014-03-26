#ifndef RefinementLine_hpp
#define RefinementLine_hpp

#include <string>
#include <vector>

#include "Geometry/ReferenceGeometry.hpp"
#include "Geometry/PhysicalGeometry.hpp"
#include "Geometry/PointPhysical.hpp"
#include "Geometry/RefinementGeometry.hpp"

namespace Geometry
{
    class RefinementLine : public Geometry::RefinementGeometry
    {
    public:
        typedef PointPhysical                PointPhysicalT;
        typedef PhysicalGeometry             PhysicalGeometryT;
        typedef ReferenceGeometry            ReferenceGeometryT;
        typedef std::vector<PointPhysicalT>     VectorOfPointPhysicalsT;
        typedef std::vector<unsigned int>       VectorOfIndicesT;

        //! Constructors.
        RefinementLine(const ReferenceGeometryT* const referenceGeometry,
                       const PhysicalGeometryT* const physicalGeometry)
                 : referenceGeometry_(referenceGeometry), physicalGeometry_(physicalGeometry)
        {
          std::cout << "RefinementLine(referenceGeometry, physicalGeometry)\n";
        }

        RefinementLine(const RefinementLine& other)
                 : referenceGeometry_(other.referenceGeometry_), 
                   physicalGeometry_(other.physicalGeometry_)
        {
          std::cout << "RefinementQuadrilateral(other)\n";
        }

        virtual ~RefinementLine()
        {}

        //! For debugging and checkpointing: a human-readable name.
        virtual std::string getName() const
        { 
          return "RefinementLine";
        }

        //---------------------- Refinement definitions -----------------------------------------
        
        //! Number of new nodes due to a refinement.
        virtual DimT nrOfNewNodes(int refineType) const;
        
        //! Get all physical nodes: existing nodes and new nodes to be added.
        virtual void getAllNodes(int refineType, VectorOfPointPhysicalsT& nodes) const;
        
        //! Number of sub-elements due to the refinement.
        virtual DimT nrOfSubElements(int refineType) const;

        //! Assembly nodes for sub-element.
        virtual void subElementLocalNodeIndices(int refineType, DimT iSubElement, VectorOfIndicesT& LocalNodeIdx) const;
        
        //! Local indices pairs of sub-elements connected by a sub-Internal Face.
        virtual void adjacentSubElementsPairs(int refineType,
                        VectorOfIndicesT& elemIdx1, VectorOfIndicesT& localFaceIdx1,
                        VectorOfIndicesT& elemIdx2, VectorOfIndicesT& localFaceIdx2) const;

        //! Number of sub-elements on a parent's face.
        virtual DimT nrOfSubElementsOnFace(int refineType, DimT faLocalIndex) const;

        //! Get sub-elements' local index on a parent's face.
        virtual void subElementsOnFace(int refineType, DimT faLocalIndex, VectorOfIndicesT& localSubElemIdx) const;
        
        //! Get sub-face's local face number of on a parent's face.
        virtual DimT getLocalSubFaceNr(int refineType, DimT localFaceNr, DimT subElementIdx) const;

    private:
        RefinementLine() : referenceGeometry_(NULL), physicalGeometry_(NULL) {}

        //! type of refinement to be applied. 

        //! The physicalGeometry object contains pointers to the actual physical nodes, and
        //! a container of global node indexes.
        const PhysicalGeometryT*        physicalGeometry_;

        //! The corresponding referenceGeometry object
        const ReferenceGeometryT* const referenceGeometry_;
    };
}

#endif
