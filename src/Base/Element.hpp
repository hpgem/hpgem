//----------------------------------------------------------------
#ifndef Element_hpp
#define Element_hpp
//----------------------------------------------------------------
#include "Base/ElementData.hpp"
#include "Base/BasisFunctionSet.hpp"
#include "Geometry/ElementGeometry.hpp"
#include "Geometry/PointPhysical.hpp"
#include "Geometry/ReferenceSquare.hpp"
#include "Geometry/PointPhysical.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRule.hpp"
#include "Base/ElementCacheData.hpp"

#include <vector>
#include <iostream>

namespace Base
{
    class Element: public Geometry::ElementGeometry,
                   public ElementData
    {
    public:
        typedef Geometry::PointPhysical                PointPhysicalT;
        typedef Geometry::PointReference               PointReferenceT;
        typedef Geometry::ReferenceGeometry            ReferenceGeometryT;
        typedef Geometry::MappingReferenceToPhysical   MappingReferenceToPhysicalT;
        typedef Geometry::ElementGeometry              ElementGeometryT;
        typedef unsigned int                                PointIndexT;
        typedef unsigned int                                UId;
        typedef std::vector<PointPhysicalT>                 VectorOfPhysicalPointsT;
        typedef std::vector<PointIndexT>                    VectorOfPointIndexesT;
        typedef Base::ElementCacheData                 CacheT;
        typedef Base::BasisFunctionSet                 BasisFunctionSetT;
        typedef QuadratureRules::GaussQuadratureRule   GaussQuadratureRuleT;
        typedef Base::ElementData                      ElementDataT;
        typedef std::vector<CacheT>                         VecCacheT;
        typedef LinearAlgebra::NumericalVector              SolutionVector;
        
    
    public:
        
        Element(const VectorOfPointIndexesT& globalNodeIndexes,
                const BasisFunctionSetT* const basisFunctionSet,
                const VectorOfPhysicalPointsT& allNodes,
                unsigned int nrOfUnkowns,
                unsigned int nrOfTimeLevels,
                unsigned int nrOfBasisFunc,
                unsigned int id,
                unsigned int numberOfElementMatrices=0,
                unsigned int numberOfElementVectors=0);
    
        Element(const Element& other);
        
        ~Element();
	

        unsigned int                    getID()const;

        unsigned int                    getID();

        void                            setQuadratureRulesWithOrder(unsigned int quadrROrder);

        void                            setGaussQuadratureRule(GaussQuadratureRuleT* const quadR);

        void							setDefaultBasisFunctionSet(BasisFunctionSetT* bFSet);

        const GaussQuadratureRuleT*     getGaussQuadratureRule() const;

        VecCacheT&                      getVecCacheData();
        
        double                          basisFunction(unsigned int i, const PointReferenceT& p) const;

		///\brief returns the value of the i-th basisfunction at point p in ret
		void                            basisFunction(unsigned int i, const PointReferenceT& p, NumericalVector& ret) const;
 
            /// jDir=0 means x, and etc.
        double                          basisFunctionDeriv(unsigned int i, unsigned int jDir, const PointReferenceT& p) const;
            //unsigned int                    getNumberOfDegreesOfFreedom()const;
            //unsigned int                    getNumberOfDegreesOfFreedom();
        
        ///\brief the all directions in one go edition of basisFunctionDeriv. Also applies the scaling gained from transforming to the reference element.
        void                           basisFunctionDeriv(unsigned int i,const PointReferenceT& p, NumericalVector& ret) const;
	
		///\brief returns the curl of the i-th basisfunction at point p in ret
		///\TODO other differential operators are not supported
		void                            basisFunctionCurl(unsigned int i, const PointReferenceT& p, NumericalVector& ret) const;
	
        void                            getSolution(unsigned int timeLevel, const PointReferenceT& p, SolutionVector& solution) const;
        
        void                            initialiseSolution(unsigned int timeLevel, unsigned int solutionId, const SolutionVector& solution);
        
    public:
            /// Output operator.
        friend ostream& operator<<(ostream& os, const Element& element)
        {
            os << '(' ;
            const Geometry::ElementGeometry& elemG = static_cast<const Geometry::ElementGeometry&>(element);
            operator<<(os, elemG);
            os<<endl;
            return os;
        }
    
    
    private:
        GaussQuadratureRuleT*                       quadratureRule_;
        const BasisFunctionSetT*                    basisFunctionSet_;
        VecCacheT                                   vecCacheData_;
        UId                                         id_;
        double                                      orderCoeff_;
    };
}

#endif
