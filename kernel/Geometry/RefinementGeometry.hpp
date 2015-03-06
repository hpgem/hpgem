/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team
 
 This code is distributed using BSD 3-Clause License. A copy of which can found below.
 
 
 Copyright (c) 2014, University of Twente
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef RefinementGeometry_hpp
#define RefinementGeometry_hpp

#include <string>
#include <vector>
             
namespace Geometry
{
	class PointPhysical;

    class RefinementGeometry
    {
    public:
        using PointPhysicalT = PointPhysical;
        using VectorOfPointPhysicalsT = std::vector<PointPhysicalT>;
        using VectorOfIndicesT = std::vector<std::size_t>;       

      public:
        virtual ~RefinementGeometry()
        {}

        //! \brief For debugging and checkpointing: a human-readable name.
        virtual std::string getName() const = 0;
        
        //---------------------- Refinement status -----------------------------------------
        
        //! \brief Get refinement type.
        std::size_t getRefineType() const { return refineType_; }

        //! \brief Set refinement type.
        void setRefineType(std::size_t refineType) { refineType_ = refineType; }

        //! \brief Unset refinement type.
        void unsetRefineType(std::size_t refineType) { refineType_ = -1; }

        //! \brief Get refinement type applied to the parent resulting this object.
        std::size_t getAppliedRefineType() const { return appliedRefineType_; }

        //! \brief Set refinement type applied to the parent resulting this object.
        void setAppliedRefineType(std::size_t appliedRefineType) { appliedRefineType_ = appliedRefineType; }

        //! \brief Is this element being refined?
        bool isBeingRefined() const { return beingRefined_; }

        //! \brief Mark: this element is being refined.
        void setBeingRefinedOn() { beingRefined_ = true; }

        //! \brief Unmark: this element is being refined.
        void setBeingRefinedOff() { beingRefined_ = false; }

        //---------------------- Refinement definitions -----------------------------------------

        //! \brief Number of new nodes due to the refinement that should be added to the vector of localNodeIndices
        virtual std::size_t nrOfNewNodes() const = 0;

        //! \brief Get all physical nodes: existing nodes and new nodes to be added.
        virtual void getAllNodes(VectorOfPointPhysicalsT& nodes) const = 0;

        //! \brief Number of sub-elements due to the refinement
        virtual std::size_t nrOfSubElements() const = 0;

        //! \brief Assembly nodes for sub-element
        virtual void subElementLocalNodeIndices(std::size_t, VectorOfIndicesT& LocalNodeIdx) const = 0;

        //! \brief Local indices pairs of sub-elements connected by a sub-Internal Face
        virtual void adjacentSubElementsPairs(
                        VectorOfIndicesT& elemIdx1, VectorOfIndicesT& localFaceIdx1,
                        VectorOfIndicesT& elemIdx2, VectorOfIndicesT& localFaceIdx2) const = 0;

        //! \brief Number of sub-elements on a particular parent's face.
        virtual std::size_t nrOfSubElementsOnFace(std::size_t faLocalIndex) const = 0;

        //! \brief Get sub-elements' local index on a particular parent's face.
        virtual void subElementsOnFace(std::size_t faLocalIndex, VectorOfIndicesT& LocalNodeIdx) const = 0;

        //! \brief Get sub-face's local face number of on a particular parent's face.
        virtual std::size_t getLocalSubFaceNr(std::size_t localFaceNr, std::size_t subElementIdx) const = 0;

      protected:
        //! \brief Default constructor.
        RefinementGeometry() 
            : refineType_(-1), 
              appliedRefineType_(-1), 
              beingRefined_(false)
        { std::cout << "RefinementGeometry()  "; }  

        //! \brief Copy constructor.
        RefinementGeometry(const RefinementGeometry& other) 
            : refineType_(other.refineType_), 
              appliedRefineType_(other.appliedRefineType_), 
              beingRefined_(other.beingRefined_)
        { std::cout << "RefinementGeometry()  "; }  

      private:
	//! this element is a result of previous refinement of this type 
        std::size_t     appliedRefineType_;
	
	//! type of refinement to be applied
        std::size_t     refineType_;
	
	//! a flag for marking unfinished refinement proccess
        bool    beingRefined_;
    };
}
#endif
