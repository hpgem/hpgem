namespace Base
{
    template<unsigned int DIM>
    ElementData<DIM>::ElementData(
        unsigned int timeLevels,
        unsigned int nrOfUnkowns,
        unsigned int nrOfBasisFunctions):
    timeLevels_(timeLevels),
    nrOfUnkowns_(nrOfUnkowns),
    nrOfBasisFunctions_(nrOfBasisFunctions)
    {
        std::vector<std::vector<double> >
        expansionCoefficients_(timeLevels_,std::vector<double>(nrOfUnkowns_*nrOfBasisFunctions_));
    }

    template<unsigned int DIM>
    std::vector<double>
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
            return expansionCoefficients_[timeLevel][unknown*basisFunction];
        }
        else
        {
            throw "Error: Asked for a time level, or unknown, greater than the amount of time levels";
        }
    }
    
    template<unsigned int DIM>
    int ElementData<DIM>::getNrOfUnknows(){return nrOfUnkowns_;}
        
    template<unsigned int DIM>
    int ElementData<DIM>::getNrOfBasisFunctions(){return nrOfBasisFunctions_;}
    
    
    
}
