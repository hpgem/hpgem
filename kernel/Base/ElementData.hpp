/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team
 
 This code is distributed using BSD 3-Clause License. A copy of which can found below.
 
 
 Copyright (c) 2014, University of Twente
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

//----------------------------------------------------------------
#ifndef ElementData_hpp
#define ElementData_hpp
//----------------------------------------------------------------
#include <vector>
#include "LinearAlgebra/Matrix.hpp"
#include "LinearAlgebra/NumericalVector.hpp"

namespace LinearAlgebra
{
  class NumericalVector;
}

struct UserElementData;

namespace Base
{

  class ElementData
  {
    /*!
     * ElementData contains the data for every element in the vector elementData_.
     * This is a two dimensional vector where the first dimensions is the time level,
     * and the second the number of unknowns times the number of basis functions.
     */
  public:
    typedef typename LinearAlgebra::Matrix               MatrixT;
    typedef typename std::vector<LinearAlgebra::Matrix>  VectorOfMatrices;
  public:

    ElementData(unsigned int timeLevels, unsigned nrOfUnkowns, unsigned int nrOfBasisFunctions, unsigned int nrOfElementMatrixes = 0, unsigned int nrOfElementVectors = 0);    
    
    virtual ~ElementData()
    {
    }

    ///Set/update the element matrix. Routines in hpGEM will assume that for every element, expansioncoefficients of unknowns belonging the the same basisfunctions
    ///appear consecutively in the matrix.
    void setElementMatrix(const LinearAlgebra::Matrix&, int matrixID=0);

    virtual void getElementMatrix(LinearAlgebra::Matrix&, int matrixID=0) const;

    ///If it not appropriate to use the timeleveldata for your vector information (for example because it is the source term in a time dependent problem)
    void setElementVector(const LinearAlgebra::NumericalVector&, int vectorID=0);

    virtual void getElementVector(LinearAlgebra::NumericalVector&, int vectorID=0) const;

    LinearAlgebra::Matrix& getTimeLevelDataMatrix(std::size_t timeLevel);

    /// Specify a time level index, return a vector containing the data for that time level.
    virtual const LinearAlgebra::NumericalVector    getTimeLevelData(size_t timeLevel) const;

    /// Specify a time level index, an unknown (as solutionId), set the data for that unknown
    void setTimeLevelData(unsigned int timeLevel, unsigned int solutionId, const LinearAlgebra::NumericalVector& unknown);
    void setTimeLevelData(unsigned int timeLevel, const LinearAlgebra::NumericalVector& unknown);

    /// Specify a time level index, an unknown and a basis function nr, return data (double)
    virtual double getData(unsigned int timeLevel, unsigned int unknown, unsigned int basisFunction) const;

    /// Specify a time level index, an unknown and a basis function nr, set the data (double)
    void setData(unsigned int timeLevel, unsigned int unknown, unsigned int basisFunction, double val);

    virtual int getNrOfUnknows() const;

    virtual int getNrOfBasisFunctions() const;

    //this needs to store information about all variables, so it needs to be a matrix (?)
    virtual const LinearAlgebra::NumericalVector& getResidue() const;

    void setResidue(LinearAlgebra::NumericalVector& residue);

    void setUserData(UserElementData* data);

    virtual UserElementData* getUserData() const
    {
      return userData_;
    }

  protected:
    void setNumberOfBasisFunctions(unsigned int number);

  private:
    size_t              timeLevels_;
    size_t              nrOfUnkowns_;
    size_t             nrOfBasisFunctions_;

    ///Stores the expansion coefficients    
    VectorOfMatrices          expansionCoefficients_;

    ///Stores the result of an element integration
    LinearAlgebra::NumericalVector           residue_;

    ///Stores polymorphic pointer to UserDefined Data, internally not used! Used only outside of the Kernel!!!
    UserElementData*          userData_;

    ///Stores element matrix(es) for this element
    VectorOfMatrices          elementMatrix_;

    std::vector<LinearAlgebra::NumericalVector> elementVector_;
  } ;
}
#endif
