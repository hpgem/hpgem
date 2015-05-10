/* 
 * File:   midPoint.cpp
 * Author: irana
 * 
 * Created on December 23, 2014, 1:43 PM
 */

#include "MidPoint.h"

namespace Base
{
    MidPoint::MidPoint()
    {
        order_ = 2;
        numStages_ = 2;
        
        //make a_
        std::vector<double> aRow;
        a_.push_back(aRow);
        aRow =
        {   0.5};
        a_.push_back(aRow);
        
        //make b_ and c_
        b_ =
        {   0.0, 1.0};
        c_ =
        {   0.0, 0.5};

    }

    std::size_t MidPoint::order() const
    {
        return order_;
    }
    
    std::size_t MidPoint::numStages() const
    {
        return numStages_;
    }
    
    double MidPoint::a(std::size_t i, std::size_t j) const
    {
        return a_[i][j];
    }
    
    double MidPoint::b(std::size_t i) const
    {
        return b_[i];
    }
    
    double MidPoint::c(std::size_t i) const
    {
        return c_[i];
    }
}
