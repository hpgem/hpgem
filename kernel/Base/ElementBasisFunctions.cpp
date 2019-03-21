/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found below.


 Copyright (c) 2019, University of Twente
 All rights reserved.

 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */


#include "ElementBasisFunctions.h"

#include "Logger.h"

namespace Base
{
    void ElementBasisFunctions::validatePositions() const
    {
        std::size_t unknowns = getNumberOfUnknowns();
        std::size_t numberOfSets = sets_->size();
        for (std::size_t i = 0; i < unknowns; ++i)
        {
            for (std::size_t j = 0; j < setPositions_[i].size(); ++j)
            {
                int position = setPositions_[i][j];
                if (position != -1) // -1 is used to signal the empty set
                {
                    logger.assert_debug(position < numberOfSets, "Invalid position");
                    logger.assert_debug(sets_->at(position) != nullptr, "Null pointer set");
                }
            }
        }
    }

    std::pair<const BasisFunctionSet*, std::size_t> ElementBasisFunctions::getBasisFunctionSetAndIndex(
            size_t index, size_t unknown) const
    {
        assertValidUnknown(unknown, true);
        if(unknown == LEGACY_BEHAVIOUR)
        {
            // Indirectly asserts constant dofs
            getNumberOfBasisFunctions();
            unknown = 0;
        }
        // Actual selection of the basis function
        std::size_t originalIndex = index;
        for(int setIndex : setPositions_[unknown])
        {
            if(setIndex != -1)
            {
                std::size_t size = sets_->at(setIndex)->size();
                if(index < size)
                {
                    return std::make_pair(sets_->at(setIndex).get(), index);
                }
                else
                {
                    index -= size;
                }
            }
        }
        logger.assert_always(false, "Asked for basis function % but there are only %", originalIndex,
                getNumberOfBasisFunctions(unknown));
        return {};
    }

    std::size_t ElementBasisFunctions::getMaximumOrder() const
    {
        std::size_t order = 0;
        for(std::size_t i = 0; i < getNumberOfUnknowns(); ++i)
        {
            for(int pos : setPositions_[i])
            {
                if(pos != -1)
                {
                    order = std::max(order, sets_->at(pos)->getOrder());
                }
            }
        }
        return order;
    }

    void ElementBasisFunctions::clearBasisFunctionPosition(std::size_t unknown)
    {
        assertValidUnknown(unknown, false);
        setPositions_[unknown].clear();
    }

    void ElementBasisFunctions::registerBasisFunctionPosition(std::size_t unknown, std::size_t place, std::size_t position)
    {
        assertValidUnknown(unknown, false);
        logger.assert_debug(position < sets_->size(), "Position % beyond the end of the set of size %",
                position, sets_->size());
        if(setPositions_[unknown].size() <= place)
        {
            setPositions_[unknown].resize(1 + place, -1);
        }
        setPositions_[unknown][place] = position;
    }
}
