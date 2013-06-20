/*
 * GlobalData.hpp
 *
 *  Created on: Feb 3, 2013
 *      Author: nicorivas
 */

#ifndef GLOBALDATA_H_
#define GLOBALDATA_H_

namespace Base
{
    struct GlobalData
    {
        virtual ~GlobalData(){}
    
        unsigned int numberOfUnknowns_;
        unsigned int numberOfTimeLevels_;
        
    };
};

#endif /* GLOBALDATA_H_ */
