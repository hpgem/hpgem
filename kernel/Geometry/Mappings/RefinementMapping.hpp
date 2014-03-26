#ifndef RefinementMapping_hpp
#define RefinementMapping_hpp

#include <iostream>
#include <string>

#include "Geometry/PointReference.hpp"
#include "LinearAlgebra/Matrix.hpp"
             
namespace Geometry
{
    class RefinementMapping
    {
    public:
        typedef unsigned int                    DimT;
        typedef PointReference             PointReferenceT;

        //! Default constructor.
        RefinementMapping() {}  

        virtual ~RefinementMapping()
        {}

        //---------------------- Refinement mappings -----------------------------------------

        //! Transform a reference point using refinement mapping
        virtual void refinementTransform(int refineType, int subElementIdx, 
                      const PointReferenceT& p, PointReferenceT& pMap) const = 0;

        //! Transformation matrix of this refinement when located on the LEFT side
        virtual void getRefinementMappingMatrixL(int refineType, int subElementIdx, 
                    LinearAlgebra::Matrix& Q) const = 0;

        //! Transformation matrix of this refinement when located on the RIGHT side
        virtual void getRefinementMappingMatrixR(int refineType, int subElementIdx, 
                    LinearAlgebra::Matrix& Q) const = 0;

        //! Refinement mapping on codim1 for a given refinement on codim0
        //! Note: this should also applied on other dimensions
        virtual void getCodim1RefinementMappingMatrixL(int refineType, DimT subElementIdx, 
                                DimT faLocalIndex, LinearAlgebra::Matrix& Q) const = 0;

        //! Refinement mapping on codim1 for a given refinement on codim0
        //! Note: this should also applied on other dimensions
        virtual void getCodim1RefinementMappingMatrixR(int refineType, DimT subElementIdx, 
                                DimT faLocalIndex, LinearAlgebra::Matrix& Q) const = 0;
    };
}
#endif
