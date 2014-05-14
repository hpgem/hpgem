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
#ifndef Face_hpp
#define Face_hpp
//----------------------------------------------------------------
#include "Geometry/GlobalNamespaceGeometry.hpp"
#include "Geometry/FaceGeometry.hpp"
#include "Base/FaceCacheData.hpp"
#include "Base/Element.hpp"
#include "Base/FaceData.hpp"

namespace Base
{
    /// Face consists of FaceGeometry and probably FaceData, if needed. FaceGeometry holds all FaceReference related data and appropriate mappings

    class Face: public Geometry::FaceGeometry,public FaceData
    {

    public:

        typedef Base::Element                                       ElementT;
        typedef Geometry::ElementGeometry                           ElementGeometryT;
        typedef typename Geometry::FaceGeometry::LocalFaceNrType    LocalFaceNrTypeT;
        typedef typename Base::FaceCacheData                        CacheT;
        typedef std::vector<CacheT>                                      VecCacheT;
        typedef Geometry::FaceGeometry                              FaceGeometryT;
        typedef QuadratureRules::GaussQuadratureRule              FaceQuadratureRule;
        
    public:
        
        Face(ElementT*  ptrElemL, const LocalFaceNrTypeT& localFaceNumL, ElementT* ptrElemRight, const LocalFaceNrTypeT& localFaceNumR,int faceID,unsigned int numberOfElementMatrixes=0,unsigned int numberOfFaceVectors=0);
        
        virtual ~Face(){}

        Face(ElementT* ptrElemL, const LocalFaceNrTypeT& localFaceNumL, const Geometry::FaceType&  ftype,int faceID, unsigned int numberOfFaceMatrixes=0, unsigned int numberOfFaceVectors=0);

        //void            setPtrElementLeft( ElementT* value);
        

        //void            setPtrElementRight( ElementT* value);

        /// Return the pointer to the left element.
        virtual ElementT*       getPtrElementLeft()     {return elementLeft_;}

        /// Return the pointer to the right element, NULL if inexistent for boundaries.
        virtual ElementT*       getPtrElementRight()    {return elementRight_;}
        
        virtual const ElementT*       getPtrElementLeft()const     {return elementLeft_;}
        
            /// Return the pointer to the right element, NULL if inexistent for boundaries.
        virtual const ElementT*       getPtrElementRight()const    {return elementRight_;}
        
        void            createQuadratureRules();
    
        void setGaussQuadratureRule(FaceQuadratureRule* quadratureRule)
        {
            quadratureRule_ = quadratureRule;
        }

        virtual const FaceQuadratureRule* getGaussQuadratureRule() const
        {
            return quadratureRule_;
        }
        
        virtual bool             isInternal()const;

        virtual VecCacheT&       getVecCacheData() { return vecCacheData_; }

        virtual double                          basisFunction(unsigned int i, const Geometry::PointReference& p) const;

		///\brief returns the value of the i-th basisfunction at point p in ret
        virtual void                            basisFunction(unsigned int i, const Geometry::PointReference& p, NumericalVector& ret) const;

        virtual void                            basisFunctionNormal(unsigned int i, const LinearAlgebra::NumericalVector& normal, const Geometry::PointReference& p, NumericalVector& ret) const;

        /// jDir=0 means x, and etc.
        virtual double                          basisFunctionDeriv(unsigned int i, unsigned int jDir, const Geometry::PointReference& p) const;

		///\brief the all directions in one go edition of basisFunctionDeriv. Also applies the scaling gained from transforming to the reference element.
        virtual void                            basisFunctionDeriv(unsigned int i,const Geometry::PointReference& p, NumericalVector& ret) const;

        virtual void                            basisFunctionCurl(unsigned int i, const Geometry::PointReference& p, NumericalVector& ret) const;

        virtual int                             getNrOfBasisFunctions() const;

        virtual int                             getLocalNrOfBasisFunctions() const{return nrOfConformingDOFOnTheFace_;}

		void                            setLocalNrOfBasisFunctions(int number){nrOfConformingDOFOnTheFace_=number;}

		virtual int getID()const{return faceID_;}

    protected:

		///\brief default constructor - for use with wrapper classes
		Face():FaceData(0,0,0),FaceGeometry(),elementLeft_(NULL),elementRight_(NULL),quadratureRule_(NULL){}

    private:
         ElementT*                                 elementLeft_;
         ElementT*                                 elementRight_;
        const FaceQuadratureRule*                             quadratureRule_;
        VecCacheT                                       vecCacheData_;

        unsigned int                                    nrOfConformingDOFOnTheFace_;
        int                                             faceID_;
    };
};
#endif
