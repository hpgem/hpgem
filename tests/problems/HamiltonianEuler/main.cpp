//
//  main.cpp
//  
//
//  Created by Shavarsh Nurijanyan on 6/21/13.
//
//

#include <iostream>
#include <cstring>

using namespace std;

#include "HEuler.hpp"

int main()
{
    const string fileName="blabla.out";
    HEulerConfigurationData config(fileName, HEulerConfigurationData::COMPRESSIBLE_PERIODIC);
    
    return 0;
}