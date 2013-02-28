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
    class GlobalData
    {
        public:
	        GlobalData();
	        virtual ~GlobalData();
        private:
	        int numberOfUnknowns;
	        int numberOfTimeLevels;
    };
};

#endif /* GLOBALDATA_H_ */
