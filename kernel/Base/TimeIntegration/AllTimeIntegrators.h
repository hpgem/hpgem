/* 
 * File:   AllTimeIntegrators.h
 * Author: irana
 *
 * Created on December 19, 2014, 11:31 AM
 */

#ifndef ALLTIMEINTEGRATORS_HPP
#define	ALLTIMEINTEGRATORS_HPP

#include "ButcherTableau.h"

namespace Base
{
    class ButcherTableau;
    
    /**
     * Storage class for all the time integrators. If you add an integrator,
     * make sure to also add it here in the constructor.
     * If you want to solve an ODE with a Butcher Tableau, this is the appropriate place 
     * to get one.
     */
    class AllTimeIntegrators
    {
    public:
        static AllTimeIntegrators& Instance();

        ButcherTableau* getRule(std::size_t order, std::size_t numStages);

        ~AllTimeIntegrators() = default;
        
    private:
        AllTimeIntegrators();

        AllTimeIntegrators(const AllTimeIntegrators& integrators) = delete;
        AllTimeIntegrators(AllTimeIntegrators&& integrators) = delete;
        
        std::vector<ButcherTableau*> vecOfIntegrators_;
        
    };
}

#endif	/* ALLTIMEINTEGRATORS_HPP */

