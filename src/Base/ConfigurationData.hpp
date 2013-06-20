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
    
        ConfigurationData(){}
        
        virtual ~ConfigurationData(){}
    
    
        unsigned int numberOfUnknowns_;
        unsigned int numberOfTimeLevels_;
        unsigned int numberOfBasisFunctions_;
    };
}

#endif /* CONFIGURATIONDATA_H_ */
