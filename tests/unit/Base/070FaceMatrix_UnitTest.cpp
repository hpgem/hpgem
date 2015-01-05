// Test the class FaceMatrix.

#include <cassert>
#include <iostream>

#include "Base/FaceMatrix.hpp"


int main()
{
    std::cout << "Test if FaceMatrix works.\n";
    
    int nDOFLeft = 3;
    int nDOFRight = 5;
    
    Base::FaceMatrix F(nDOFLeft, nDOFRight);
    for(int i=0; i < nDOFLeft+nDOFRight; i++)
    {
        for(int j=0; j < nDOFLeft+nDOFRight; j++)
        {
            F(i,j) = i*10+j;
        }
    }
    
    Base::FaceMatrix F2(F);
    F2.axpy(100,F);
    
    Base::FaceMatrix F3(F);
    F3 *= 100;
    F3 += F;
    
    Base::Side iS=Base::Side::LEFT;
    Base::Side jS=Base::Side::LEFT;
    int iVB=0;
    int jVB=0;
    int z=0;
    for(int i=0; i < nDOFLeft+nDOFRight; i++)
    {
        if(i < nDOFLeft)
        {
            iS=Base::Side::LEFT;
            iVB=i;
        }
        else
        {
            iS=Base::Side::RIGHT;
            iVB=i-nDOFLeft;
        }
        
        for(int j=0; j < nDOFLeft+nDOFRight; j++)
        {
            if(j < nDOFLeft)
            {
                jS=Base::Side::LEFT;
                jVB=j;
            }
            else
            {
                jS=Base::Side::RIGHT;
                jVB=j-nDOFLeft;
            }
            
            z=i*10+j;
            z=z*100+z;
            
            assert(F2.getEntireMatrix()(i,j) == z);
            assert(F2.getElementMatrix(iS,jS)(iVB,jVB) == z);
            assert(F3(i,j) == z);
        }
    }
    
    //std::cout << "Matrix F2:\n" << F2.getEntireMatrix() << "\n";
    return 0;
}