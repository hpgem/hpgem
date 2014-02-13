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
        
        Face(ElementT*  ptrElemL, const LocalFaceNrTypeT& localFaceNumL, ElementT* ptrElemRight, const LocalFaceNrTypeT& localFaceNumR,unsigned int numberOfElementMatrixes=0,unsigned int numberOfFaceVectors=0);
        
        virtual ~Face(){}

        Face(ElementT* ptrElemL, const LocalFaceNrTypeT& localFaceNumL, const Geometry::FaceType&  ftype, unsigned int numberOfFaceMatrixes=0, unsigned int numberOfFaceVectors=0);

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
#endif
