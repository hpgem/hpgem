#include "AllTimeIntegrators.h"
#include "ButcherTableau.h"
#include "RK4Methods.h"
#include "ForwardEuler.h"
#include "MidPoint.h"

namespace Base
{
    AllTimeIntegrators::AllTimeIntegrators()
    {
        vecOfIntegrators_.push_back(&ForwardEuler::instance());
        vecOfIntegrators_.push_back(&MidPoint::instance());
        vecOfIntegrators_.push_back(&RK4_4::instance());
        //Should we sort the rules, or is it okay like this?
    }
    
    AllTimeIntegrators& AllTimeIntegrators::Instance()
    {
        static AllTimeIntegrators theInstance;
        return theInstance;
    }
    
    ButcherTableau* AllTimeIntegrators::getRule(std::size_t order, std::size_t numStages)
    {
        for (ButcherTableau* rule : vecOfIntegrators_)
        {
            if (rule->order() == order && rule->numStages() == numStages)
            {
                return rule;
            }
        }
        throw "Could not find the Runge Kutta method you're looking for.";
    }
}
