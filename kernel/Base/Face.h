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

//----------------------------------------------------------------
#ifndef Face_h
#define Face_h
//----------------------------------------------------------------
#include "Base/FaceData.h"
#include "Base/FaceMatrix.h"
#include "Geometry/FaceGeometry.h"
#include "Base/Element.h"
#include "L2Norm.h"
#include "TreeEntry.h"

namespace QuadratureRules
{
    class GaussQuadratureRule;
}

namespace Base
{
    class Element;
    struct FaceCacheData;
    
    /// Face consists of FaceGeometry and probably FaceData, if needed. 
    ///FaceGeometry holds all FaceReference related data and appropriate mappings
    //class is final as a reminder that the virtual default destructor should be added once something inherits from this class
    class Face final: public Geometry::FaceGeometry, public FaceData
    {
    public:
        
        using FaceQuadratureRule = QuadratureRules::GaussQuadratureRule;

        Face(Element* ptrElemL, const std::size_t& localFaceNumL,
                Element* ptrElemRight, const std::size_t& localFaceNumR,
                std::size_t faceID, std::size_t numberOfFaceMatrixes = 0, 
                std::size_t numberOfFaceVectors = 0);
        
        Face(Element* ptrElemL, const std::size_t& localFaceNumL,
        const Geometry::FaceType& ftype, std::size_t faceID, 
        std::size_t numberOfFaceMatrixes = 0, std::size_t numberOfFaceVectors = 0);
                
        ///copy constructor should not be used: if adjacent elements are the same,
        ///then the Face already exists and there is no need for another, if 
        ///adjacent elements are different, the copy is not really a copy
        Face(const Face &other) = delete;
        Face& operator=(const Face &other) = delete;
        
        ///Copy constructor with new elements. It makes a copy of the face, but 
        ///with new elements assigned to it.
        Face(const Face& other, Element* elementL, const std::size_t localFaceL, Element* elementR, const std::size_t localFaceR);
        
        /// Return the pointer to the left element.
        Element* getPtrElementLeft()
        {
            return elementLeft_;
        }
        
        /// Return the pointer to the right element, nullptr if inexistent for boundaries.
        Element* getPtrElementRight()
        {
            return elementRight_;
        }
        
        const Element* getPtrElementLeft() const
        {
            return elementLeft_;
        }
        
        /// Return the pointer to the right element, nullptr if inexistent for boundaries.
        const Element* getPtrElementRight() const
        {
            return elementRight_;
        }
        
        /// \brief Return the pointer to the element on side iSide.
        const Element* getPtrElement(Side iSide) const
        {
            if (iSide == Side::LEFT)
            {
                return elementLeft_;
            }
            else
            {
                return elementRight_;
            }
        }
        
        /// \brief Return the pointer to the element on the opposite side of iSide.
        const Element* getPtrElementOpposite(Side iSide) const
        {
            if (iSide == Side::LEFT)
            {
                return elementRight_;
            }
            else
            {
                return elementLeft_;
            }
        }

        /// \brief Get a list of pointers to the adjacent nodes
        const std::vector<Base::Node *> getNodesList() const;
        
        ///get the root element that is the (indirect) parent of one of the adjacent elements
        Element* getRootElement();

        /// \brief Create a quadrature for this face based on the quadrature rules of adjacent elements.
        void createQuadratureRules();

        void setGaussQuadratureRule(FaceQuadratureRule* quadratureRule)
        {
            logger.assert(quadratureRule!=nullptr, "Invalid quadrature rule passed");
            quadratureRule_ = quadratureRule;
        }
        
        /// \brief Get a pointer to the quadrature rule used to do integration on this face.
        FaceQuadratureRule* getGaussQuadratureRule() const
        {
            return quadratureRule_;
        }
        
        /// \brief Get the value of the basis function (corresponding to index i) at the physical point corresponding to reference point p.
        template<std::size_t DIM>
        double basisFunction(std::size_t i, const Geometry::PointReference<DIM>& p) const;
        template<std::size_t DIM>
        double basisFunction(std::size_t i, const Geometry::PointReference<DIM>& p, std::size_t unknown) const;
        
        ///\brief returns the value of the i-th basis function at point p in ret
        template<std::size_t DIM>
        void basisFunction(std::size_t i, const Geometry::PointReference<DIM>& p, LinearAlgebra::SmallVector<DIM + 1>& ret) const;
        template<std::size_t DIM>
        void basisFunction(std::size_t i, const Geometry::PointReference<DIM>& p, LinearAlgebra::SmallVector<DIM + 1>& ret, std::size_t unknown) const;
        
        /// \brief Returns the value of the basis function (corresponding to element function index iBasisFunction) on the element at side iSide at the physical point corresponding to reference point p.
        template<std::size_t DIM>
        double basisFunction(Side iSide, std::size_t iBasisFunction, const Geometry::PointReference<DIM>& p) const;
        template<std::size_t DIM>
        double basisFunction(Side iSide, std::size_t iBasisFunction, const Geometry::PointReference<DIM>& p, std::size_t unknown) const;
        
        /// \brief Returns the value of the i-th basis function multiplied by the outward pointing normal vector. The value is computed at the physical point corresponding to reference point p.
        template<std::size_t DIM>
        LinearAlgebra::SmallVector<DIM + 1> basisFunctionNormal(std::size_t i, const LinearAlgebra::SmallVector<DIM + 1>& normal, const Geometry::PointReference<DIM>& p) const;
        template<std::size_t DIM>
        LinearAlgebra::SmallVector<DIM + 1> basisFunctionNormal(std::size_t i, const LinearAlgebra::SmallVector<DIM + 1>& normal, const Geometry::PointReference<DIM>& p, std::size_t unknown) const;

        /// \brief Returns the physical normal vector multiplied by the basis function (corresponding to element function index iBasisFunction) on the element at side iSide. The value is computed at the physical point corresponding to reference point p.
        template<std::size_t DIM>
        LinearAlgebra::SmallVector<DIM + 1> basisFunctionNormal(Side iSide, std::size_t iBasisFunction, const LinearAlgebra::SmallVector<DIM + 1>& normal, const Geometry::PointReference<DIM>& p) const;
        template<std::size_t DIM>
        LinearAlgebra::SmallVector<DIM + 1> basisFunctionNormal(Side iSide, std::size_t iBasisFunction, const LinearAlgebra::SmallVector<DIM + 1>& normal, const Geometry::PointReference<DIM>& p, std::size_t unknown) const;
        
        /// \brief Returns the (physical) derivative in direction jDir of the physical basis function (corresponding to index i) at the physical point corresponding to reference point p.
        template<std::size_t DIM>
        double basisFunctionDeriv(std::size_t i, std::size_t jDir, const Geometry::PointReference<DIM>& p) const;
        template<std::size_t DIM>
        double basisFunctionDeriv(std::size_t i, std::size_t jDir, const Geometry::PointReference<DIM>& p, std::size_t unknown) const;

        ///\brief The "all directions in one go"-edition of basisFunctionDeriv. Warning: contrary to some old implementations, this will NOT apply the scaling due to the transformation to a reference element
        template<std::size_t DIM>
        LinearAlgebra::SmallVector<DIM + 1> basisFunctionDeriv(std::size_t i, const Geometry::PointReference<DIM>& p) const;
        template<std::size_t DIM>
        LinearAlgebra::SmallVector<DIM + 1> basisFunctionDeriv(std::size_t i, const Geometry::PointReference<DIM>& p, std::size_t unknown) const;

        /// \brief Returns the (physical) gradient of the physical basis function (corresponding to element function index iBasisFunction) on the element at side iSide. The gradient is computed at the physical point corresponding to reference point p.
        template<std::size_t DIM>
        LinearAlgebra::SmallVector<DIM + 1> basisFunctionDeriv(Side iSide, std::size_t iBasisFunction, const Geometry::PointReference<DIM>& p) const;
        template<std::size_t DIM>
        LinearAlgebra::SmallVector<DIM + 1> basisFunctionDeriv(Side iSide, std::size_t iBasisFunction, const Geometry::PointReference<DIM>& p, std::size_t unknown) const;

        /// \brief Returns the (physical) curl of the physical basis function (corresponding to index i) at the physical point corresponding to reference point p.
        template<std::size_t DIM>
        LinearAlgebra::SmallVector<DIM + 1> basisFunctionCurl(std::size_t i, const Geometry::PointReference<DIM>& p) const;
        template<std::size_t DIM>
        LinearAlgebra::SmallVector<DIM + 1> basisFunctionCurl(std::size_t i, const Geometry::PointReference<DIM>& p, std::size_t unknown) const;

        /// \brief Returns the sum of the number of basisfunctions of the adjacent elements.
        std::size_t getNumberOfBasisFunctions() const;
        std::size_t getNumberOfBasisFunctions(std::size_t unknown) const;

        std::size_t getLocalNumberOfBasisFunctions() const
        {
            std::size_t number = numberOfConformingDOFOnTheFace_[0];
            for(std::size_t index : numberOfConformingDOFOnTheFace_)
                logger.assert(index == number, "number of basis functions is different for different unknown");
            return numberOfConformingDOFOnTheFace_[0];
        }
        
        std::size_t getLocalNumberOfBasisFunctions(std::size_t unknown) const
        {
            logger.assert(unknown < numberOfConformingDOFOnTheFace_.size(),
                          "Asking for unknown % but there are only %", unknown, numberOfConformingDOFOnTheFace_.size());
            return numberOfConformingDOFOnTheFace_[unknown];
        }

        std::size_t getTotalLocalNumberOfBasisFunctions() const
        {
            std::size_t result = 0;
            for (auto nbasis : numberOfConformingDOFOnTheFace_)
                result += nbasis;
            return result;
        }
        
        void setLocalNumberOfBasisFunctions(std::size_t number, std::size_t unknown = 0)
        {
            logger.assert(unknown < numberOfConformingDOFOnTheFace_.size(),
                          "Setting for unknown % but there are only %", unknown, numberOfConformingDOFOnTheFace_.size());
            numberOfConformingDOFOnTheFace_[unknown] = number;
        }
        
        ///\deprecated Does not conform naming conventions, use getNumberOfBasisFunctions instead.
        std::size_t getNrOfBasisFunctions() const
        {
            return getNumberOfBasisFunctions();
        }

        ///\deprecated Does not conform naming conventions, use getLocalNumberOfBasisFunctions instead.
        std::size_t getLocalNrOfBasisFunctions() const
        {
            return getLocalNumberOfBasisFunctions();
        }
        
        ///\deprecated Does not conform naming conventions, use setLocalNumberOfBasisFunctions instead.
        void setLocalNrOfBasisFunctions(std::size_t number)
        {
            setLocalNumberOfBasisFunctions(number);
        }
        
        std::size_t getID() const
        {
            return faceID_;
        }
        
        /// Specify a time integration vector id, return a vector containing the data for that time integration vector.
        LinearAlgebra::MiddleSizeVector getTimeIntegrationVector(std::size_t timeIntegrationVectorId, std::size_t unknown = 0) const;

        /// \brief Convert the side of the face, the index corresponding to the scalar basis function (scalarBasisFunctionId) and the index corresponding to the variable (varId) to a single index (faceBasisFunctionId).
        std::size_t convertToSingleIndex(Side side, std::size_t scalarBasisFunctionId, std::size_t varId = 0) const;
        
        /// \brief Convert the index of the basis (vector)function of the face (faceBasisFunctionId) to the corresponding side of the face.
        Side getSide(std::size_t faceBasisFunctionId) const;
        
        /// \brief Convert the index of the basis (vector)function of the face (faceBasisFunctionId) to the index of the corresponding element basis (vector)function (elementBasisFunctionId).
        std::size_t getElementBasisFunctionId(std::size_t faceBasisFunctionId) const;

        void addElement(Element* ptrElementR, std::size_t localFaceNumberR);

        void setPositionInTree(const TreeEntry<Face*>* position) {
            logger.assert(position->getData() == this, "Trying to set the position of another face as this face");
            positionInTheTree_ = position;
        }

        const TreeEntry<Face*>* getPositionInTree() const {
            return positionInTheTree_;
        }
    private:
        

        Element* elementLeft_;
        Element* elementRight_;
        FaceQuadratureRule* quadratureRule_;

        const TreeEntry<Face*>* positionInTheTree_;

        std::vector<std::size_t> numberOfConformingDOFOnTheFace_;
        std::size_t faceID_;
    };
}
#include "FaceCacheData.h"
namespace Base
{
    template<std::size_t DIM>
    double Face::basisFunction(std::size_t i, const Geometry::PointReference<DIM>& p) const
    {
        logger.assert(i<getNumberOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", i, getNumberOfBasisFunctions());
        std::size_t numberOfBasisFuncs = getPtrElementLeft()->getNumberOfBasisFunctions();
        if (i < numberOfBasisFuncs)
        {
            return getPtrElementLeft()->basisFunction(i, mapRefFaceToRefElemL(p));
        }
        else
        {
            return getPtrElementRight()->basisFunction(i - numberOfBasisFuncs, mapRefFaceToRefElemR(p));
        }
    }
    
    template<std::size_t DIM>
    double Face::basisFunction(std::size_t i, const Geometry::PointReference<DIM>& p, std::size_t unknown) const
    {
        logger.assert(i<getNumberOfBasisFunctions(unknown), "Asked for basis function %, but there are only % basis functions", i, getNumberOfBasisFunctions(unknown));
        std::size_t numberOfBasisFuncs = getPtrElementLeft()->getNumberOfBasisFunctions(unknown);
        if (i < numberOfBasisFuncs)
        {
            return getPtrElementLeft()->basisFunction(i, mapRefFaceToRefElemL(p), unknown);
        }
        else
        {
            return getPtrElementRight()->basisFunction(i - numberOfBasisFuncs, mapRefFaceToRefElemR(p), unknown);
        }
    }

    template<std::size_t DIM>
    void Face::basisFunction(std::size_t i, const Geometry::PointReference<DIM>& p, LinearAlgebra::SmallVector<DIM + 1>& ret) const
    {
        logger.assert(i<getNumberOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", i, getNumberOfBasisFunctions());
        std::size_t n(getPtrElementLeft()->getNumberOfBasisFunctions());
        if (i < n)
        {
            getPtrElementLeft()->basisFunction(i, mapRefFaceToRefElemL(p), ret);
        }
        else
        {
            getPtrElementRight()->basisFunction(i - n, mapRefFaceToRefElemR(p), ret);
        }
    }
    
    template<std::size_t DIM>
    void Face::basisFunction(std::size_t i, const Geometry::PointReference<DIM>& p, LinearAlgebra::SmallVector<DIM + 1>& ret, std::size_t unknown) const
    {
        logger.assert(i<getNumberOfBasisFunctions(unknown), "Asked for basis function %, but there are only % basis functions", i, getNumberOfBasisFunctions(unknown));
        std::size_t n(getPtrElementLeft()->getNumberOfBasisFunctions(unknown));
        if (i < n)
        {
            getPtrElementLeft()->basisFunction(i, mapRefFaceToRefElemL(p), ret, unknown);
        }
        else
        {
            getPtrElementRight()->basisFunction(i - n, mapRefFaceToRefElemR(p), ret, unknown);
        }
    }

    /// \param[in] iSide The index corresponding to the side of the face.
    /// \param[in] iBasisFunction The index corresponding to the basis function.
    /// \param[in] p The reference point on the reference element.
    template<std::size_t DIM>
    double Face::basisFunction(Side iSide, std::size_t iBasisFunction, const Geometry::PointReference<DIM>& p) const
    {
        if (iSide == Side::LEFT)
        {
            logger.assert(iBasisFunction < getPtrElementLeft()->getNumberOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", iBasisFunction, getPtrElementLeft()->getNumberOfBasisFunctions());
            return getPtrElementLeft()->basisFunction(iBasisFunction, mapRefFaceToRefElemL(p));
        }
        else
        {
            logger.assert(isInternal(), "boundary faces only have a \"left\" element");
            logger.assert(iBasisFunction < getPtrElementRight()->getNumberOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", iBasisFunction, getPtrElementRight()->getNumberOfBasisFunctions());
            return getPtrElementRight()->basisFunction(iBasisFunction, mapRefFaceToRefElemR(p));
        }
    }
    
    /// \param[in] iSide The index corresponding to the side of the face.
    /// \param[in] iBasisFunction The index corresponding to the basis function.
    /// \param[in] p The reference point on the reference element.
    template<std::size_t DIM>
    double Face::basisFunction(Side iSide, std::size_t iBasisFunction, const Geometry::PointReference<DIM>& p, std::size_t unknown) const
    {
        if (iSide == Side::LEFT)
        {
            logger.assert(iBasisFunction < getPtrElementLeft()->getNumberOfBasisFunctions(unknown), "Asked for basis function %, but there are only % basis functions", iBasisFunction, getPtrElementLeft()->getNumberOfBasisFunctions(unknown));
            return getPtrElementLeft()->basisFunction(iBasisFunction, mapRefFaceToRefElemL(p), unknown);
        }
        else
        {
            logger.assert(isInternal(), "boundary faces only have a \"left\" element");
            logger.assert(iBasisFunction < getPtrElementRight()->getNumberOfBasisFunctions(unknown), "Asked for basis function %, but there are only % basis functions", iBasisFunction, getPtrElementRight()->getNumberOfBasisFunctions(unknown));
            return getPtrElementRight()->basisFunction(iBasisFunction, mapRefFaceToRefElemR(p), unknown);
        }
    }

    template<std::size_t DIM>
    LinearAlgebra::SmallVector<DIM + 1> Face::basisFunctionNormal(std::size_t i, const LinearAlgebra::SmallVector<DIM + 1>& normal, const Geometry::PointReference<DIM>& p) const
    {
        logger.assert(i<getNumberOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", i, getNumberOfBasisFunctions());
        LinearAlgebra::SmallVector<DIM + 1> ret;
        std::size_t n = getPtrElementLeft()->getNumberOfBasisFunctions();
        if (i < n)
        {
            ret = normal;
            ret *= getPtrElementLeft()->basisFunction(i, mapRefFaceToRefElemL(p)) / Base::L2Norm(normal);
        }
        else
        {
            ret = normal;
            ret *= -getPtrElementRight()->basisFunction(i - n, mapRefFaceToRefElemR(p)) / Base::L2Norm(normal);
        }
        return ret;
    }
    
    template<std::size_t DIM>
    LinearAlgebra::SmallVector<DIM + 1> Face::basisFunctionNormal(std::size_t i, const LinearAlgebra::SmallVector<DIM + 1>& normal, const Geometry::PointReference<DIM>& p, std::size_t unknown) const
    {
        logger.assert(i<getNumberOfBasisFunctions(unknown), "Asked for basis function %, but there are only % basis functions", i, getNumberOfBasisFunctions(unknown));
        LinearAlgebra::SmallVector<DIM + 1> ret;
        std::size_t n = getPtrElementLeft()->getNumberOfBasisFunctions(unknown);
        if (i < n)
        {
            ret = normal;
            ret *= getPtrElementLeft()->basisFunction(i, mapRefFaceToRefElemL(p), unknown) / Base::L2Norm(normal);
        }
        else
        {
            ret = normal;
            ret *= -getPtrElementRight()->basisFunction(i - n, mapRefFaceToRefElemR(p), unknown) / Base::L2Norm(normal);
        }
        return ret;
    }

    /// \param[in] iSide The index corresponding to the side of the face.
    /// \param[in] iBasisFunction The index corresponding to the basis function.
    /// \param[in] normal The normal vector (pointing outwards with respect to the element on the left side).
    /// \param[in] p The reference point on the reference element.
    template<std::size_t DIM>
    LinearAlgebra::SmallVector<DIM + 1> Face::basisFunctionNormal(Side iSide, std::size_t iBasisFunction, const LinearAlgebra::SmallVector<DIM + 1>& normal, const Geometry::PointReference<DIM>& p) const
    {
        if (iSide == Side::LEFT)
        {
            logger.assert(iBasisFunction < getPtrElementLeft()->getNumberOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", iBasisFunction, getPtrElementLeft()->getNumberOfBasisFunctions());
            return getPtrElementLeft()->basisFunction(iBasisFunction, mapRefFaceToRefElemL(p)) * normal / Base::L2Norm(normal);
        }
        else
        {
            logger.assert(isInternal(), "boundary faces only have a \"left\" element");
            logger.assert(iBasisFunction < getPtrElementRight()->getNumberOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", iBasisFunction, getPtrElementRight()->getNumberOfBasisFunctions());
            return -getPtrElementRight()->basisFunction(iBasisFunction, mapRefFaceToRefElemR(p)) * normal / Base::L2Norm(normal);
        }
    }
    
    /// \param[in] iSide The index corresponding to the side of the face.
    /// \param[in] iBasisFunction The index corresponding to the basis function.
    /// \param[in] normal The normal vector (pointing outwards with respect to the element on the left side).
    /// \param[in] p The reference point on the reference element.
    template<std::size_t DIM>
    LinearAlgebra::SmallVector<DIM + 1> Face::basisFunctionNormal(Side iSide, std::size_t iBasisFunction, const LinearAlgebra::SmallVector<DIM + 1>& normal, const Geometry::PointReference<DIM>& p, std::size_t unknown) const
    {
        if (iSide == Side::LEFT)
        {
            logger.assert(iBasisFunction < getPtrElementLeft()->getNumberOfBasisFunctions(unknown), "Asked for basis function %, but there are only % basis functions", iBasisFunction, getPtrElementLeft()->getNumberOfBasisFunctions(unknown));
            return getPtrElementLeft()->basisFunction(iBasisFunction, mapRefFaceToRefElemL(p), unknown) * normal / Base::L2Norm(normal);
        }
        else
        {
            logger.assert(isInternal(), "boundary faces only have a \"left\" element");
            logger.assert(iBasisFunction < getPtrElementRight()->getNumberOfBasisFunctions(unknown), "Asked for basis function %, but there are only % basis functions", iBasisFunction, getPtrElementRight()->getNumberOfBasisFunctions(unknown));
            return -getPtrElementRight()->basisFunction(iBasisFunction, mapRefFaceToRefElemR(p), unknown) * normal / Base::L2Norm(normal);
        }
    }

    template<std::size_t DIM>
    double Face::basisFunctionDeriv(std::size_t i, std::size_t jDir, const Geometry::PointReference<DIM>& p) const
    {
        logger.assert(i<getNumberOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", i, getNumberOfBasisFunctions());
        std::size_t n = getPtrElementLeft()->getNumberOfBasisFunctions();
        if (i < n)
        {
            return getPtrElementLeft()->basisFunctionDeriv(i, jDir, mapRefFaceToRefElemL(p));
        }
        else
        {
            return getPtrElementRight()->basisFunctionDeriv(i - n, jDir, mapRefFaceToRefElemR(p));
        }
    }
    
    template<std::size_t DIM>
    double Face::basisFunctionDeriv(std::size_t i, std::size_t jDir, const Geometry::PointReference<DIM>& p, std::size_t unknown) const
    {
        logger.assert(i<getNumberOfBasisFunctions(unknown), "Asked for basis function %, but there are only % basis functions", i, getNumberOfBasisFunctions(unknown));
        std::size_t n = getPtrElementLeft()->getNumberOfBasisFunctions(unknown);
        if (i < n)
        {
            return getPtrElementLeft()->basisFunctionDeriv(i, jDir, mapRefFaceToRefElemL(p), unknown);
        }
        else
        {
            return getPtrElementRight()->basisFunctionDeriv(i - n, jDir, mapRefFaceToRefElemR(p), unknown);
        }
    }

    template<std::size_t DIM>
    LinearAlgebra::SmallVector<DIM + 1> Face::basisFunctionDeriv(std::size_t i, const Geometry::PointReference<DIM>& p) const
    {
        logger.assert(i<getNumberOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", i, getNumberOfBasisFunctions());
        std::size_t n = getPtrElementLeft()->getNumberOfBasisFunctions();
        if (i < n)
        {
            return getPtrElementLeft()->basisFunctionDeriv(i, mapRefFaceToRefElemL(p));
        }
        else
        {
            return getPtrElementRight()->basisFunctionDeriv(i - n, mapRefFaceToRefElemR(p));
        }
    }
    
    template<std::size_t DIM>
    LinearAlgebra::SmallVector<DIM + 1> Face::basisFunctionDeriv(std::size_t i, const Geometry::PointReference<DIM>& p, std::size_t unknown) const
    {
        logger.assert(i<getNumberOfBasisFunctions(unknown), "Asked for basis function %, but there are only % basis functions", i, getNumberOfBasisFunctions(unknown));
        std::size_t n = getPtrElementLeft()->getNumberOfBasisFunctions(unknown);
        if (i < n)
        {
            return getPtrElementLeft()->basisFunctionDeriv(i, mapRefFaceToRefElemL(p), unknown);
        }
        else
        {
            return getPtrElementRight()->basisFunctionDeriv(i - n, mapRefFaceToRefElemR(p), unknown);
        }
    }

    /// \param[in] iSide The index corresponding to the side of the face.
    /// \param[in] iBasisFunction The index corresponding to the basis function.
    /// \param[in] p The reference point on the reference element.
    template<std::size_t DIM>
    LinearAlgebra::SmallVector<DIM + 1> Face::basisFunctionDeriv(Side iSide, std::size_t iBasisFunction, const Geometry::PointReference<DIM>& p) const
    {
        if (iSide == Side::LEFT)
        {
            logger.assert(iBasisFunction < getPtrElementLeft()->getNumberOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", iBasisFunction, getPtrElementLeft()->getNumberOfBasisFunctions());
            return getPtrElementLeft()->basisFunctionDeriv(iBasisFunction, mapRefFaceToRefElemL(p));
        }
        else
        {
            logger.assert(isInternal(), "boundary faces only have a \"left\" element");
            logger.assert(iBasisFunction < getPtrElementRight()->getNumberOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", iBasisFunction, getPtrElementRight()->getNumberOfBasisFunctions());
            return getPtrElementRight()->basisFunctionDeriv(iBasisFunction, mapRefFaceToRefElemR(p));
        }
    }
    
    /// \param[in] iSide The index corresponding to the side of the face.
    /// \param[in] iBasisFunction The index corresponding to the basis function.
    /// \param[in] p The reference point on the reference element.
    template<std::size_t DIM>
    LinearAlgebra::SmallVector<DIM + 1> Face::basisFunctionDeriv(Side iSide, std::size_t iBasisFunction, const Geometry::PointReference<DIM>& p, std::size_t unknown) const
    {
        if (iSide == Side::LEFT)
        {
            logger.assert(iBasisFunction < getPtrElementLeft()->getNumberOfBasisFunctions(unknown), "Asked for basis function %, but there are only % basis functions", iBasisFunction, getPtrElementLeft()->getNumberOfBasisFunctions(unknown));
            return getPtrElementLeft()->basisFunctionDeriv(iBasisFunction, mapRefFaceToRefElemL(p), unknown);
        }
        else
        {
            logger.assert(isInternal(), "boundary faces only have a \"left\" element");
            logger.assert(iBasisFunction < getPtrElementRight()->getNumberOfBasisFunctions(unknown), "Asked for basis function %, but there are only % basis functions", iBasisFunction, getPtrElementRight()->getNumberOfBasisFunctions(unknown));
            return getPtrElementRight()->basisFunctionDeriv(iBasisFunction, mapRefFaceToRefElemR(p), unknown);
        }
    }

    template<std::size_t DIM>
    LinearAlgebra::SmallVector<DIM + 1> Face::basisFunctionCurl(std::size_t i, const Geometry::PointReference<DIM>& p) const
    {
        logger.assert(i<getNumberOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", i, getNumberOfBasisFunctions());
        std::size_t numberOfBasisFunctionsLeft = getPtrElementLeft()->getNumberOfBasisFunctions();
        if (i < numberOfBasisFunctionsLeft)
        {
            return getPtrElementLeft()->basisFunctionCurl(i, mapRefFaceToRefElemL(p));
        }
        else
        {
            return getPtrElementRight()->basisFunctionCurl(i - numberOfBasisFunctionsLeft, mapRefFaceToRefElemR(p));
        }
    }
    
    template<std::size_t DIM>
    LinearAlgebra::SmallVector<DIM + 1> Face::basisFunctionCurl(std::size_t i, const Geometry::PointReference<DIM>& p, std::size_t unknown) const
    {
        logger.assert(i<getNumberOfBasisFunctions(unknown), "Asked for basis function %, but there are only % basis functions", i, getNumberOfBasisFunctions(unknown));
        std::size_t numberOfBasisFunctionsLeft = getPtrElementLeft()->getNumberOfBasisFunctions(unknown);
        if (i < numberOfBasisFunctionsLeft)
        {
            return getPtrElementLeft()->basisFunctionCurl(i, mapRefFaceToRefElemL(p), unknown);
        }
        else
        {
            return getPtrElementRight()->basisFunctionCurl(i - numberOfBasisFunctionsLeft, mapRefFaceToRefElemR(p), unknown);
        }
    }
}
#endif
