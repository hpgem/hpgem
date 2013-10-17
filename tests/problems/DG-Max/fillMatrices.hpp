#ifndef FILLMATRICES_HPP
#define FILLMATRICES_HPP

template <unsigned int DIM>
class hpGemUIExtentions;

/**
 * this class is a wrapper for the global assembly routines, so a choise about using the IP-method or the Brezzi method can be made on a high level
 */
class matrixFiller{
public:
    /// the actual global assembly routine
    virtual void fillMatrixes(hpGemUIExtentions<3>* matrixContainer)=0;
};

/**
 * this class provides a global assembly structure that is suitable for use with the IP-method
 */
class matrixFillerIP:public matrixFiller{
public:
    virtual void fillMatrixes(hpGemUIExtentions<3>* matrixContainer);
};

/**
 * this class provides a global assembly stucture that is suitable for the Brezzi-flux
 */
class matrixFillerBR:public matrixFiller{
public:
    virtual void fillMatrixes(hpGemUIExtentions<3>* matrixContainer);
};

#endif