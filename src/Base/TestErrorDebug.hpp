//-------------------------------------------------------------------
// MT Julianto, Sat Feb 16 20:10:07 WET 2013
// Trap error on DEBUG mode
//-------------------------------------------------------------------
#ifndef TestErrorDebug_hpp
#define TestErrorDebug_hpp

#include <cassert>

#ifdef DEBUG
#define TestErrorDebug(boolCondition, errMessage) assert((errMessage,boolCondition));
#else
#ifndef NDEBUG
#define TestErrorDebug(boolCondition, errMessage) assert((errMessage,boolCondition));
#else
#define TestErrorDebug(boolCondition, errMessage)
#endif
#endif

#endif  // TestErrorDebug_hpp
