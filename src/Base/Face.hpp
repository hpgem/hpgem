//----------------------------------------------------------------
#ifndef Face_hpp
#define Face_hpp
//----------------------------------------------------------------
#include "Geometry/GlobalNamespaceGeometry.hpp"
#include "Geometry/FaceGeometry.hpp"
#include "Base/FaceCacheData.hpp"
#include "Base/Element.hpp"

namespace Base
{
    /// Face consists of FaceGeometry and probably FaceData, if needed. FaceGeometry holds all FaceReference related data and appropriate mappings

    template <unsigned int DIM>
    class Face: public Geometry::FaceGeometry<DIM>
    {

    public:

        typedef Base::Element<DIM>                                       ElementT;
        typedef Geometry::ElementGeometry<DIM>                           ElementGeometryT;
        typedef typename Geometry::FaceGeometry<DIM>::LocalFaceNrType    LocalFaceNrTypeT;
        typedef typename Base::FaceCacheData<DIM>                        CacheT;
        typedef std::vector<CacheT>                                      VecCacheT;
        typedef Geometry::FaceGeometry<DIM>                              FaceGeometryT;
        typedef QuadratureRules::GaussQuadratureRule<DIM-1>              FaceQuadratureRule;
        
    public:
        
        Face(ElementT*  ptrElemL, const LocalFaceNrTypeT& localFaceNumL, ElementT* ptrElemRight, const LocalFaceNrTypeT& localFaceNumR);
        
        virtual ~Face(){}

        Face(ElementT* ptrElemL, const LocalFaceNrTypeT& localFaceNumL, const Geometry::FaceType&  ftype);

        void            setPtrElementLeft(const ElementT* value);
        

        void            setPtrElementRight(const ElementT* value);

        /// Return the pointer to the left element.
        const ElementT*       getPtrElementLeft()     {return elementLeft_;}

        /// Return the pointer to the right element, NULL if inexistent for boundaries.
        const ElementT*       getPtrElementRight()    {return elementRight_;}
        
        const ElementT*       getPtrElementLeft()const     {return elementLeft_;}
        
            /// Return the pointer to the right element, NULL if inexistent for boundaries.
        const ElementT*       getPtrElementRight()const    {return elementRight_;}
        
        void            createQuadratureRules();
    
        void setGaussQuadratureRule(FaceQuadratureRule* quadratureRule)
        {
            quadratureRule_ = quadratureRule;
        }

        FaceQuadratureRule* getGaussQuadratureRule() const
        {
            return quadratureRule_;
        }
        
        bool             isInternal()const;

        VecCacheT&       getVecCacheData() { return vecCacheData_; }

    private:
        const ElementT*                                 elementLeft_;
        const ElementT*                                 elementRight_;
        FaceQuadratureRule*                             quadratureRule_;
        VecCacheT                                       vecCacheData_;
    };
};
#include "Face_Impl.hpp"
#endif
