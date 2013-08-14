/*
 * ConfigurationData.hpp
 *
 *  Created on: Feb 3, 2013
 *      Author: nicorivas
 */

#ifndef CONFIGURATIONDATA_H_
#define CONFIGURATIONDATA_H_

namespace Base
{
    struct ConfigurationData
    {
    
        ConfigurationData(unsigned int numberOfUnknowns, unsigned int numberOfBasisFunctions, unsigned int  numberOfTimeLevels=1):
            numberOfUnknowns_(numberOfUnknowns),
            numberOfBasisFunctions_(numberOfBasisFunctions),
            numberOfTimeLevels_(numberOfTimeLevels)
        {}
        
        virtual ~ConfigurationData(){}
    
    
        unsigned int numberOfUnknowns_;
        unsigned int numberOfBasisFunctions_;
        unsigned int numberOfTimeLevels_;
    };
}

#endif /* CONFIGURATIONDATA_H_ */
