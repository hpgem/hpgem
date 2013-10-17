
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
        expansionCoefficients_(timeLevels_),
        userData_(NULL)
    {
        for (typename VectorOfMatrices::iterator cit=expansionCoefficients_.begin(); cit!=expansionCoefficients_.end(); ++cit)
            cit->resize(nrOfUnkowns_, nrOfBasisFunctions);
    }

    template<unsigned int DIM>
    const LinearAlgebra::Matrix&
    ElementData<DIM>::getTimeLevelData(unsigned int timeLevel) const
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
    ElementData<DIM>::getData(unsigned int timeLevel, unsigned int unknown, unsigned int basisFunction) const
    {
        if (timeLevel < timeLevels_ && unknown < nrOfUnkowns_ * nrOfBasisFunctions_)
        {
            return expansionCoefficients_[timeLevel](unknown, basisFunction);
        }
        else
        {
            throw "Error: Asked for a time level, or unknown, greater than the amount of time levels";
        }
    }
    template<unsigned int DIM>
    void
    ElementData<DIM>::setData(unsigned int timeLevel, unsigned int unknown, unsigned int basisFunction, double val)
    {
        if (timeLevel < timeLevels_ && unknown < nrOfUnkowns_ * nrOfBasisFunctions_)
        {
            expansionCoefficients_[timeLevel](unknown, basisFunction)=val;
        }
        else
        {
            throw "Error: Asked for a time level, or unknown, greater than the amount of time levels";
        }
    }
    

    
        ///Rewrite with swap!!! and for all variables immediately
    template<unsigned int DIM>
    void
    ElementData<DIM>::setTimeLevelData(unsigned int timeLevel, unsigned int solutionId, const LinearAlgebra::NumericalVector& unknown)
    {
        if (timeLevel < timeLevels_ && solutionId < nrOfUnkowns_)
        {
                //cout << "came here with "<< "solutionId="<<solutionId<< ",unknown="<<unknown[0]<<endl;
            LinearAlgebra::Matrix& mat = expansionCoefficients_[timeLevel];
                // cout << mat<<"before setting"<<endl;
    
            for (int i = 0; i < unknown.size(); ++i)
            {
                mat(solutionId, i) = unknown[i];
//                cout << "was here"<<solutionId<<endl;
            }
            
                // cout << mat<<"after setting"<<endl;
        }
        else
        {
            throw "Error: Asked for a time level, or unknown, greater than the amount of time levels";
        }
    }
    
    template<unsigned int DIM>
    int
    ElementData<DIM>::getNrOfUnknows() const
    {
        return nrOfUnkowns_;
    }
        
    template<unsigned int DIM>
    int
    ElementData<DIM>::getNrOfBasisFunctions() const
    {
        return nrOfBasisFunctions_;
    }
    
    template<unsigned int DIM>
    typename ElementData<DIM>::VectorOfDoubles&
    ElementData<DIM>::getResidue() const
    {
        return residue_;
    }
    
    template<unsigned int DIM>
    void
    ElementData<DIM>::setResidue(VectorOfDoubles& residue)
    {
        residue_=residue;
    }
    template<unsigned int DIM>
    void
    ElementData<DIM>::setUserData(UserElementData* data)
    {
        userData_=data;
    }
}
