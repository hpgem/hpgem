//----------------------------------------------------------------
#ifndef FaceData_hpp
#define FaceData_hpp
//----------------------------------------------------------------
#include <vector>
#include "FaceCacheData.hpp"

namespace Base
{
    template <unsigned int DIM>
    class FaceData
    {
      public:
        typedef FaceCacheData<DIM> CacheT;
        typedef std::vector<CacheT> VecCacheT;

      public:
        FaceData() {}
        
        VecCacheT& getVecCacheData() const
        {
          return vecCacheData_;
        }

        virtual ~FaceData() {;}

      private:
        VecCacheT vecCacheData_;
    };
};
#include "FaceData_Impl.hpp"
#endif