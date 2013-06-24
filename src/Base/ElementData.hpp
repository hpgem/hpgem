//----------------------------------------------------------------
#ifndef ElementData_hpp
#define ElementData_hpp
//----------------------------------------------------------------
#include <vector>
#include "Base/Element.hpp"
#include "Base/UserElementData.hpp"
#include "LinearAlgebra/Matrix.hpp"

namespace Base
{
    template <unsigned int DIM>
	class ElementData
	{
        /*!
         * ElementData contains the data for every element in the vector elementData_.
         * This is a two dimensional vector where the first dimensions is the time level,
         * and the second the number of unknowns times the number of basis functions.
         */
    public:
        typedef typename std::vector<double>                 VectorOfDoubles;
        typedef typename std::vector<LinearAlgebra::Matrix>  VectorOfMatrices;
        
    public:
// MTJ: default constructor, just for a moment; to be removed later!  Soon!
        ElementData() {}
	      
        ElementData(unsigned int timeLevels, unsigned nrOfUnkowns, unsigned int nrOfBasisFunctions);

       
        virtual ~ElementData() {delete userData_;}
        
        
            /// Specify a time level index, return a vector containing the data for that time level.
        const LinearAlgebra::Matrix&    getTimeLevelData(unsigned int timeLevel) const;
        
            /// Specify a time level index, an unknown and a basis function nr, return data (double)
        double                    getData(unsigned int timeLevel, unsigned int unknown, unsigned int basisFunction) const;

        
        int                       getNrOfUnknows() const;
        
        int                       getNrOfBasisFunctions() const;
        
        VectorOfDoubles&          getResidue() const;
        
        void                      setResidue(VectorOfDoubles& residue);
        
        void                      setUserData(UserElementData* data);
        
        UserElementData*          getUserData();
            
    private:
        unsigned int              timeLevels_;
        unsigned int              nrOfUnkowns_;
        unsigned int              nrOfBasisFunctions_;
        ///Stores the expansion coefficients    
        VectorOfMatrices          expansionCoefficients_;
        ///Stores the result of an element integration
        VectorOfDoubles           residue_;
        ///Stores polymorphic pointer to UserDefined Data, internally not used! Used only outside of the Kernel!!!
        UserElementData*          userData_;
	};
}
#include "ElementData_Impl.hpp"
#endif
