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

#ifndef REFERENCEGEOMETRY_HH
#define REFERENCEGEOMETRY_HH

#include "Geometry/PointReference.hpp"
#include "Geometry/Mappings/MappingCodimensions.hpp"
#include "Geometry/Mappings/RefinementMapping.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRule.hpp"
#include <map>
#include <unordered_map>
#include "Base/BaseBasisFunction.hpp"

#include <iostream>
#include <vector>
using std::cout;
using std::endl;

template<>
class std::__1::hash<Geometry::PointReference>{
public:
	size_t operator()(const Geometry::PointReference& point) const{
		static std::__1::hash<double> hasher;
		size_t ret=0;
		for(int i=0;i<point.size();++i){
			ret^=hasher(point[i])+0x9e3779b9+(ret<<6)+(ret>>2);
		}
		return ret;
	}
};

namespace QuadratureRules
{
    // forward declaration
    class GaussQuadratureRule;
}

namespace Geometry
{
    
    enum TypeOfReferenceGeometry
    {
        POINT,
        LINE,
        TRIANGLE,
        SQUARE,
        TETRAHEDRON,
        PYRAMID,
        CUBE,
        TRIANGULARPRISM,
        HYPERCUBE
    };
    
    class ReferenceGeometry :
            public RefinementMapping,
            public MappingCodimensions
    {
    /*! \class ReferenceGeometry
     * \brief ReferenceGeometry stores a the information of a unitary geometry where the integration is made.
     * \details
     * ReferenceGeometry stores the information of the corresponding reference object, which
     * are used for integration routines. It is a pure virtual class; an interface for every
     * particular reference geometry. The information is a container of the actual reference
     * points.
     *
     * Constructors are protected, to avoid the creation of two identical physical geometries.
     * Only one of every needed type is necessary.
     */
    public:
        /// \bug this is a workaround for a g++ bug. Should read using typenames;
        typedef std::string                                         String;
        typedef unsigned int                                        IndexT;
        typedef typename Geometry::PointReference              PointReferenceT;
        typedef typename std::vector<PointReferenceT >              VectorOfReferencePointsT;
        typedef typename VectorOfReferencePointsT::iterator         iterator;
        typedef typename VectorOfReferencePointsT::const_iterator   const_iterator;
        typedef std::vector<IndexT>                                 ListOfIndexesT;

    public:

        virtual ~ReferenceGeometry(){};

        /// \brief Check whether a given point is within the ReferenceGeometry.
        virtual bool            isInternalPoint(const PointReferenceT& point) const = 0;

        /// \brief Each reference geometry knows its center of mass.
        virtual void            getCenter(PointReferenceT& point) const = 0;

        /// \brief Return number of nodes of the reference shape.
        virtual unsigned int    getNumberOfNodes() const {return points_.size();}
        TypeOfReferenceGeometry getGeometryType() const  {return geometryType_;}

        /// \brief Given a local index, return (assign to point) the corresponding node.
        virtual void            getNode(const IndexT& localIndex, PointReferenceT& node) const
                                {node = points_[localIndex];}

        virtual int             getLocalNodeIndex(int face, int node)const = 0;

        /// \brief For debugging and checkpointing: a human-readable name.
        virtual std::string     getName() const = 0;
        
        // ================================== Quadrature rules =====================================
        
        /// \brief Get a valid quadrature for this geometry.
        virtual const QuadratureRules::GaussQuadratureRule* const getGaussQuadratureRule(int order) const;

        ///\bug getBasisFunctionValue and getBasisFunctionDerivative have functionality that is completely independent from the rest of ReferenceGeometry
        ///\bug getBasisFunctionValue does some lazy initialization, so it can't be const, unless you consider the state to
        /// contain the values of all basisFunctions at all reference points
        double getBasisFunctionValue(const Base::BaseBasisFunction* function, const PointReference& p);

        double getBasisFunctionValue(const Base::BaseBasisFunction* function, const PointReference& p) const
        {return const_cast<ReferenceGeometry*>(this)->getBasisFunctionValue(function,p);}

        ///\bug getBasisFunctionDerivative does some lazy initialization, so it can't be const, unless you consider the state to
        /// contain the values of all basisFunctions at all reference points
        void getBasisFunctionDerivative(const Base::BaseBasisFunction* function, const PointReference& p,NumericalVector& ret);

        void getBasisFunctionDerivative(const Base::BaseBasisFunction* function, const PointReference& p,NumericalVector& ret) const
        {const_cast<ReferenceGeometry*>(this)->getBasisFunctionDerivative(function,p,ret);}

    protected:
        ReferenceGeometry(const TypeOfReferenceGeometry& geoT);
        ReferenceGeometry(unsigned int numberOfNodes, unsigned int DIM, const TypeOfReferenceGeometry& geoT);
        ReferenceGeometry(const ReferenceGeometry& other);
        
    protected:
        /// Container of the actual points (no reference).
        VectorOfReferencePointsT        points_;
        /// An identifier of the type of referenceGeometry, that some say shouldn't be used.
        const TypeOfReferenceGeometry   geometryType_;
        
    private:

        std::map<const Base::BaseBasisFunction*,std::unordered_map<Geometry::PointReference,double> > basisfunctionValues_;
        std::map<const Base::BaseBasisFunction*,std::unordered_map<Geometry::PointReference,NumericalVector> > basisfunctionDerivatives_;

    };
};
#endif
