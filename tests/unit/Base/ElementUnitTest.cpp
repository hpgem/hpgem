#include "Base/Element.hpp"
#include "Geometry/PointPhysical.hpp"
#include "Geometry/PhysicalGeometry.hpp"
#include "Geometry/PhysicalGeometry.hpp"
#include "Base/BasisFunctionSet.hpp"
#include "Base/AssembleBasisFunctionSet.hpp"

#include <vector>

using namespace std;

typedef std::vector<Geometry::PointPhysical<2> >        VectorOfPhysicalPointsT;
typedef Base::BasisFunctionSet<2>                       BasisFunctionSetT;
typedef std::vector<BasisFunctionSetT*>                 CollectionOfBasisFunctionSets;

int main()
{
/*
    std::vector<unsigned int> PointsIndexes(4);
    for (int i=0;i<4;i++) {PointsIndexes[i]=i;}
    
    VectorOfPhysicalPointsT points;
    
    Geometry::PointPhysical<2> point;
    
    
    //End points of a square with coordinates (0,0), (1,0), (1,1) and (0,1)
    point[0]=0;
    point[1]=0;
    points[0]=point;
    
    point[0]=1;
    points[1]=point;
    
    point[0]=1;
    points[2]=point;
    
    point[0]=0;
    points[3]=point;

    unsigned int numOfUnknowns      = 10;
    unsigned int numOfTimeLevels    = 10;
    unsigned int counter            = 10;
    
    
    CollectionOfBasisFunctionSets   collBasisFSet_;
    
    Base::BasisFunctionSet<2>* bFset1 = new Base::BasisFunctionSet<2>(2);
    Base::AssembleBasisFunctionSet_2D_Ord2_A1(*bFset1);
    collBasisFSet_.push_back(bFset1);
    
    const BasisFunctionSetT* const bf= collBasisFSet_[0];
    
    
    Base::Element<2>  myElement(PointsIndexes, bf, points, numOfUnknowns, numOfTimeLevels, 11, counter);

 */
    return 0;
}
