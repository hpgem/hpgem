/* 
 * File:   RK4Methods.h
 * Author: irana
 *
 * Created on December 19, 2014, 10:48 AM
 */

#ifndef RK4METHODS_HPP
#define	RK4METHODS_HPP

#include <vector>
#include "ButcherTableau.h"

namespace Base
{
    class RK4_4 : public ButcherTableau
    {
    public:
        static RK4_4& instance()
        {
            static RK4_4 theInstance;
            return theInstance;
        }
        
        std::size_t order() const;
        std::size_t numStages() const;
        double a(std::size_t i, std::size_t j) const;
        double b(std::size_t i) const;
        double c(std::size_t i) const;

    private:
        RK4_4();
        virtual ~RK4_4()
        {
        }
        
        std::size_t order_;
        std::size_t numStages_;
        std::vector<std::vector<double> > a_;
        std::vector<double> b_;
        std::vector<double> c_;
    };
}

#endif	/* RK4METHODS_HPP */

