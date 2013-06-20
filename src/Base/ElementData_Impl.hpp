
namespace Base
{
    template<unsigned int DIM>
    class ElementData;
    
    template<unsigned int DIM>
    ElementData<DIM>::ElementData(unsigned int timeLevels,
                                  unsigned int nrOfUnkowns,
                                  unsigned int nrOfBasisFunctions):
        timeLevels_(timeLevels),
        nrOfUnkowns_(nrOfUnkowns),
        nrOfBasisFunctions_(nrOfBasisFunctions),
        expansionCoefficients_(timeLevels_)
    {
        for (typename VectorOfMatrices::iterator cit=expansionCoefficients_.begin(); cit!=expansionCoefficients_.end(); ++cit)
            cit->resize(nrOfBasisFunctions, nrOfUnkowns_);
    }

    template<unsigned int DIM>
    LinearAlgebra::Matrix&
    ElementData<DIM>::getTimeLevelData(unsigned int timeLevel)
    {
        if (timeLevel < timeLevels_)
        {
            return expansionCoefficients_[timeLevel];
        }
        else
        {
            throw "Error: Asked for a time level greater than the amount of time levels";
        }
    }


    template<unsigned int DIM>
    double
    ElementData<DIM>::getData(unsigned int timeLevel, unsigned int unknown, unsigned int basisFunction)
    {
        if (timeLevel < timeLevels_ && unknown < nrOfUnkowns_ * nrOfBasisFunctions_)
        {
            return expansionCoefficients_[timeLevel](unknown,basisFunction);
        }
        else
        {
            throw "Error: Asked for a time level, or unknown, greater than the amount of time levels";
        }
    }
    
    template<unsigned int DIM>
    int
    ElementData<DIM>::getNrOfUnknows()
    {
        return nrOfUnkowns_;
    }
        
    template<unsigned int DIM>
    int
    ElementData<DIM>::getNrOfBasisFunctions()
    {
        return nrOfBasisFunctions_;
    }
    
    template<unsigned int DIM>
    typename ElementData<DIM>::VectorOfDoubles&
    ElementData<DIM>::getResidue()
    {
        return residue_;
    }
    
    template<unsigned int DIM>
    void
    ElementData<DIM>::setResidue(VectorOfDoubles& residue)
    {
        residue_=residue;
    }
}
