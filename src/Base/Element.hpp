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
    template <unsigned int DIM>
    class Element: public Geometry::ElementGeometry<DIM>,
                   public ElementData<DIM>
    {
    public:
        typedef Geometry::PointPhysical<DIM>                PointPhysicalT;
        typedef Geometry::PointReference<DIM>               PointReferenceT;
        typedef Geometry::ReferenceGeometry<DIM>            ReferenceGeometryT;
        typedef Geometry::MappingReferenceToPhysical<DIM>   MappingReferenceToPhysicalT;
        typedef Geometry::ElementGeometry<DIM>              ElementGeometryT;
        typedef unsigned int                                PointIndexT;
        typedef unsigned int                                UId;
        typedef std::vector<PointPhysicalT>                 VectorOfPhysicalPointsT;
        typedef std::vector<PointIndexT>                    VectorOfPointIndexesT;
        typedef Base::ElementCacheData<DIM>                 CacheT;
        typedef Base::BasisFunctionSet<DIM>                 BasisFunctionSetT;
        typedef QuadratureRules::GaussQuadratureRule<DIM>   GaussQuadratureRuleT;
        typedef Base::ElementData<DIM>                      ElementDataT;
        typedef std::vector<CacheT>                         VecCacheT;
        typedef LinearAlgebra::NumericalVector              SolutionVector;
        
    
    public:
        
        Element(const VectorOfPointIndexesT& globalNodeIndexes,
                const BasisFunctionSetT* const basisFunctionSet,
                const VectorOfPhysicalPointsT& allNodes,
                unsigned int nrOfUnkowns,
                unsigned int nrOfTimeLevels,
                unsigned int nrOfBasisFunc,
                unsigned int id);
    
        Element(const Element& other);
        
        ~Element();

        unsigned int                    getID()const;

        unsigned int                    getID();

        void                            setQuadratureRulesWithOrder(unsigned int quadrROrder);

        void                            setGaussQuadratureRule(GaussQuadratureRuleT* const quadR);

        const GaussQuadratureRuleT*     getGaussQuadratureRule() const;

        VecCacheT&                      getVecCacheData();
        
        double                          basisFunction(unsigned int i, const PointReferenceT& p) const;
 
            /// jDir=0 means x, and etc.
        double                          basisFunctionDeriv(unsigned int i, unsigned int jDir, const PointReferenceT& p) const;
            //unsigned int                    getNumberOfDegreesOfFreedom()const;
            //unsigned int                    getNumberOfDegreesOfFreedom();
        
        void                            getSolution(unsigned int timeLevel, const PointReferenceT& p, SolutionVector& solution) const;
        
        void                            initialiseSolution(unsigned int timeLevel, unsigned int solutionId, const SolutionVector& solution);
        
    public:
            /// Output operator.
        friend ostream& operator<<(ostream& os, const Element& element)
        {
            os << '(' ;
            const Geometry::ElementGeometry<DIM>& elemG = static_cast<const Geometry::ElementGeometry<DIM>&>(element);
            operator<<(os, elemG);
            os<<endl;
            return os;
        }
    
    
    private:
        const BasisFunctionSetT* const              basisFunctionSet_;
        GaussQuadratureRuleT*                       quadratureRule_;
        VecCacheT                                   vecCacheData_;
        UId                                         id_;
        double                                      orderCoeff_;
    };
}
#include "Element_Impl.hpp"

#endif
