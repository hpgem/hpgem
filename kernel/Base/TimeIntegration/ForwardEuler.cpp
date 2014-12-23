/* 
 * File:   forwardEuler.cpp
 * Author: irana
 * 
 * Created on December 23, 2014, 1:31 PM
 */

#include "ForwardEuler.hpp"

namespace Base
{
    ForwardEuler::ForwardEuler()
    {
        order_ = 1;
        numStages_ = 1;
        
        //make a_
        std::vector<double> aRow;
        a_.push_back (aRow);
        
        //make b_ and c_
        b_ = {1.0};
        c_ = {0.0};

    }
    
    std::size_t ForwardEuler::order() const
    {
        return order_;
    }

    std::size_t ForwardEuler::numStages() const
    {
        return numStages_;
    }
    
    double ForwardEuler::a(std::size_t i, std::size_t j) const
    {
        return a_[i][j];
    }
    
    double ForwardEuler::b(std::size_t i) const
    {
        return b_[i];
    }
    
    double ForwardEuler::c(std::size_t i) const
    {
        return c_[i];
    }
}
