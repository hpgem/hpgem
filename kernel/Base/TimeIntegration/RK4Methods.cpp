#include "RK4Methods.h"

namespace Base
{
    RK4_4::RK4_4()
    {
        order_ = 4;
        numStages_ = 4;
        
        //make a_
        std::vector<double> aRow;
        a_.push_back(aRow);
        aRow =
        {   0.5};
        a_.push_back(aRow);
        aRow =
        {   0.0, 0.5};
        a_.push_back(aRow);
        aRow =
        {   0.0, 0.0, 1.0};
        a_.push_back(aRow);
        
        //make b_ and c_
        b_ =
        {   1.0/6, 1.0/3, 1.0/3, 1.0/6};
        c_ =
        {   0.0, 0.5, 0.5, 1.0};

    }

    std::size_t RK4_4::order() const
    {
        return order_;
    }
    
    std::size_t RK4_4::numStages() const
    {
        return numStages_;
    }
    
    double RK4_4::a(std::size_t i, std::size_t j) const
    {
        return a_[i][j];
    }
    
    double RK4_4::b(std::size_t i) const
    {
        return b_[i];
    }
    
    double RK4_4::c(std::size_t i) const
    {
        return c_[i];
    }
}
