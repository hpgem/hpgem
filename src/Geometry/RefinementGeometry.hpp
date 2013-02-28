#ifndef RefinementGeometry_hpp
#define RefinementGeometry_hpp

#include <string>
#include <vector>

#include "Geometry/PointPhysical.hpp"
             
namespace Geometry
{
    template <unsigned int DIM>
    class RefinementGeometry
    {
    public:
        typedef unsigned int                    DimT;
        typedef PointPhysical<DIM>              PointPhysicalT;
        typedef std::vector<PointPhysicalT>     VectorOfPointPhysicalsT;
        typedef std::vector<unsigned int>       VectorOfIndicesT;

      public:
        virtual ~RefinementGeometry()
        {}

        //! \brief For debugging and checkpointing: a human-readable name.
        virtual std::string getName() const = 0;
        
        //---------------------- Refinement status -----------------------------------------
        
        //! \brief Get refinement type.
        int getRefineType() const { return refineType_; }

        //! \brief Set refinement type.
        void setRefineType(int refineType) { refineType_ = refineType; }

        //! \brief Unset refinement type.
        void unsetRefineType(int refineType) { refineType_ = -1; }

        //! \brief Get refinement type applied to the parent resulting this object.
        int getAppliedRefineType() const { return appliedRefineType_; }

        //! \brief Set refinement type applied to the parent resulting this object.
        void setAppliedRefineType(int appliedRefineType) { appliedRefineType_ = appliedRefineType; }

        //! \brief Is this element being refined?
        bool isBeingRefined() const { return beingRefined_; }

        //! \brief Mark: this element is being refined.
        void setBeingRefinedOn() { beingRefined_ = true; }

        //! \brief Unmark: this element is being refined.
        void setBeingRefinedOff() { beingRefined_ = false; }

        //---------------------- Refinement definitions -----------------------------------------

        //! \brief Number of new nodes due to the refinement that should be added to the vector of localNodeIndices
        virtual DimT nrOfNewNodes() const = 0;

        //! \brief Get all physical nodes: existing nodes and new nodes to be added.
        virtual void getAllNodes(VectorOfPointPhysicalsT& nodes) const = 0;

        //! \brief Number of sub-elements due to the refinement
        virtual DimT nrOfSubElements() const = 0;

        //! \brief Assembly nodes for sub-element
        virtual void subElementLocalNodeIndices(DimT, VectorOfIndicesT& LocalNodeIdx) const = 0;

        //! \brief Local indices pairs of sub-elements connected by a sub-Internal Face
        virtual void adjacentSubElementsPairs(
                        VectorOfIndicesT& elemIdx1, VectorOfIndicesT& localFaceIdx1,
                        VectorOfIndicesT& elemIdx2, VectorOfIndicesT& localFaceIdx2) const = 0;

        //! \brief Number of sub-elements on a parent's face.
        virtual DimT nrOfSubElementsOnFace(DimT faLocalIndex) const = 0;

        //! \brief Get sub-elements' local index on a parent's face.
        virtual void subElementsOnFace(DimT faLocalIndex, VectorOfIndicesT& LocalNodeIdx) const = 0;

        //! \brief Get sub-face's local face number of on a parent's face.
        virtual DimT getLocalSubFaceNr(DimT localFaceNr, DimT subElementIdx) const = 0;

      protected:
        //! \brief Default constructor.
        RefinementGeometry() : refineType_(-1), appliedRefineType_(-1), beingRefined_(false)
        { std::cout << "RefinementGeometry()  "; }  

      private:
        int     refineType_;
        int     appliedRefineType_;
        bool    beingRefined_;
    };
}
#endif
