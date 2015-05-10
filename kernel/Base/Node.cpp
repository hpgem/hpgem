/* 
 * File:   Node.cpp
 * Author: brinkf
 * 
 * Created on November 18, 2014, 3:13 PM
 */

#include "Node.h"
#include "Element.h"
#include "LinearAlgebra/NumericalVector.h"
#include "ElementCacheData.h"

void Base::Node::addElement(Element* element, std::size_t localNodeNr)
{
    elements_.push_back(element);
    localNodeNrs_.push_back(localNodeNr);
    element->setNode(localNodeNr, this);
}

Base::Element* Base::Node::getElement(std::size_t i)
{
    return elements_[i];
}

const Base::Element* Base::Node::getElement(std::size_t i) const
{
    return elements_[i];
}

std::size_t Base::Node::getNrOfElements() const
{
    return elements_.size();
}
