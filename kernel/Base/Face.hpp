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
         ElementT*       getPtrElementLeft()     {return elementLeft_;}

        /// Return the pointer to the right element, NULL if inexistent for boundaries.
         ElementT*       getPtrElementRight()    {return elementRight_;}
        
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

        double                          basisFunction(unsigned int i, const Geometry::PointReference& p) const;

		///\brief returns the value of the i-th basisfunction at point p in ret
		void                            basisFunction(unsigned int i, const Geometry::PointReference& p, NumericalVector& ret) const;

		void                            basisFunctionNormal(unsigned int i, const LinearAlgebra::NumericalVector& normal, const Geometry::PointReference& p, NumericalVector& ret) const;

        /// jDir=0 means x, and etc.
		double                          basisFunctionDeriv(unsigned int i, unsigned int jDir, const Geometry::PointReference& p) const;

		///\brief the all directions in one go edition of basisFunctionDeriv. Also applies the scaling gained from transforming to the reference element.
		void                            basisFunctionDeriv(unsigned int i,const Geometry::PointReference& p, NumericalVector& ret) const;

		void                            basisFunctionCurl(unsigned int i, const Geometry::PointReference& p, NumericalVector& ret) const;

		int                             getNrOfBasisFunctions() const;

		int                             getLocalNrOfBasisFunctions() const{return nrOfConformingDOFOnTheFace_;}

		void                            setLocalNrOfBasisFunctions(int number){nrOfConformingDOFOnTheFace_=number;}

		int getID()const{return faceID_;}

    private:
         ElementT*                                 elementLeft_;
         ElementT*                                 elementRight_;
        FaceQuadratureRule*                             quadratureRule_;
        VecCacheT                                       vecCacheData_;

        unsigned int                                    nrOfConformingDOFOnTheFace_;
        int                                             faceID_;
    };
};
#endif
