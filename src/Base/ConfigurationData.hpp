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
    
        ConfigurationData(unsigned int DIMension, unsigned int numberOfUnknowns, unsigned int polynomialOrder, unsigned int  numberOfTimeLevels=1):
            numberOfUnknowns_(numberOfUnknowns),
            polynomialOrder_(polynomialOrder),
            numberOfTimeLevels_(numberOfTimeLevels),
            dimension_(DIMension)
        {}
        
        virtual ~ConfigurationData(){}
    
        unsigned int dimension_;
        unsigned int numberOfUnknowns_;
        unsigned int numberOfBasisFunctions_;
        unsigned int numberOfTimeLevels_;
        unsigned int polynomialOrder_;
    };
}

#endif /* CONFIGURATIONDATA_H_ */
