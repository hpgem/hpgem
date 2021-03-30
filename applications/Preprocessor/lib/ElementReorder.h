/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2014, University of Twente
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.

 3. Neither the name of the copyright holder nor the names of its contributors
 may be used to endorse or promote products derived from this software without
 specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef HPGEM_APP_ELEMENTREORDER_H
#define HPGEM_APP_ELEMENTREORDER_H

#include <vector>
#include <string>


namespace Preprocessor {

class ElementReorder {

   public:
    // order is a vector of indices of the nodes in the mesh format in the order
    // of hpgem e.g. hpgem ordering is 0,1,2 and your format is 1,0,2
    // order={1,0,2}
    void addElementType(size_t dimension, const std::string& name,
                        const std::vector<size_t>& order);

    void reorderToHpGem(size_t dimension, std::vector<size_t>& indeces) const;

    void reorderFromHpGem(size_t dimension, std::vector<size_t>& indeces) const;

   private:
    struct Element {

        Element(size_t dimension, const std::string& name,
                const std::vector<size_t>& order)
            : dimension_(dimension), name_(name), order_(order) {
            checkOrder(order);
        }
        size_t dimension_;
        std::string name_;
        std::vector<size_t> order_;

        bool operator==(const Element& other) {
            return dimension_ == other.dimension_ && name_ == other.name_ &&
                   order_ == other.order_;
        }

       private:
        void checkOrder(std::vector<size_t> order) const;
    };

    const Element& FindElement(size_t dimension, size_t indices_size) const;

    std::vector<Element> order_per_element_;  // there are like 10 different
                                              // types so a linear search should
                                              // be fast enough
};

}  // namespace Preprocessor

#endif
