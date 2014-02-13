//----------------------------------------------------------------
#ifndef FaceData_hpp
#define FaceData_hpp
//----------------------------------------------------------------
#include <vector>
#include "FaceCacheData.hpp"
#include "LinearAlgebra/Matrix.hpp"

namespace Base
{
    class FaceData
    {
      public:
        typedef FaceCacheData CacheT;
        typedef std::vector<CacheT> VecCacheT;

      public:
        FaceData(unsigned int numberOfDOF, unsigned int numberOfFaceMatrices=0, unsigned int numberOfFaceVactors=0);
        
        void
        setFaceMatrix(const LinearAlgebra::Matrix& matrix, unsigned int matrixID=0);

        void
        getFaceMatrix(LinearAlgebra::Matrix& matrix, unsigned int matrixID=0) const;

        void
        setFaceVector(const LinearAlgebra::NumericalVector& vector, unsigned int vectorID=0);

        void
        getFaceVector(LinearAlgebra::NumericalVector& vector, unsigned int vectorID=0) const;

        const VecCacheT& getVecCacheData() const
        {
          return vecCacheData_;
        }

        virtual ~FaceData() {;}

      private:
        VecCacheT vecCacheData_;
        std::vector<LinearAlgebra::Matrix> faceMatrix_;
        std::vector<LinearAlgebra::NumericalVector> faceVector_;
    };
};
#endif
