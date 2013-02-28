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
        std::vector<std::vector<NumberT> >
        elementData_(timeLevels_,std::vector<NumberT>(nrOfUnkowns_*nrOfBasisFunctions_));
    }

    template<unsigned int DIM>
    std::vector<typename ElementData<DIM>::NumberT>
    ElementData<DIM>::getTimeLevelData(unsigned int timeLevel)
    {
        if (timeLevel < timeLevels_)
        {
            return elementData_[timeLevel];
        }
        else
        {
            throw "Error: Asked for a time level greater than the amount of time levels";
        }
    }

    template<unsigned int DIM>
    double
    ElementData<DIM>::getData(unsigned int timeLevel, unsigned int index)
    {
        if (timeLevel < timeLevels_ && index < nrOfUnkowns_ * nrOfBasisFunctions_)
        {
            return elementData_[timeLevel][index];
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
            return elementData_[timeLevel][unknown*basisFunction];
        }
        else
        {
            throw "Error: Asked for a time level, or unknown, greater than the amount of time levels";
        }
    }
}
