#include "ElementInfos.hpp"

void InvertAndTranspose(Geometry::Jacobian& orig, Geometry::Jacobian& inverse){
    //direct computation using the definitions of the inverse and the transpose
    inverse(0,0)=orig(1,1)*orig(2,2)-orig(1,2)*orig(2,1);
    inverse(1,0)=orig(2,1)*orig(0,2)-orig(2,2)*orig(0,1);
    inverse(2,0)=orig(0,1)*orig(1,2)-orig(0,2)*orig(1,1);
    inverse(0,1)=orig(1,2)*orig(2,0)-orig(1,0)*orig(2,2);
    inverse(1,1)=orig(2,2)*orig(0,0)-orig(2,0)*orig(0,2);
    inverse(2,1)=orig(0,2)*orig(1,0)-orig(0,0)*orig(1,2);
    inverse(0,2)=orig(1,0)*orig(2,1)-orig(1,1)*orig(2,0);
    inverse(1,2)=orig(2,0)*orig(0,1)-orig(2,1)*orig(0,0);
    inverse(2,2)=orig(0,0)*orig(1,1)-orig(0,1)*orig(1,0);
    inverse/=orig.determinant();
}

void FunctionCache::getFunctionValuesVector(const Base::Element* element, const PointElementReferenceT& point, std::vector< NumericalVector >& values){
    values=valueCache_[point];
    if(values.empty()){
	NumericalVector value(3);
	for(int j=0;j<element->getNrOfBasisFunctions();++j){
	    element->basisFunction(j,point,value);
	    values.push_back(value);
	} 
	valueCache_[point]=values;
    }
}

void FunctionCache::getFunctionCurlsVector(const Base::Element* element, const PointElementReferenceT& point, std::vector< NumericalVector >& curls){
    curls=curlCache_[point];
    if(curls.empty()){
	NumericalVector curl(3);
	for(int j=0;j<element->getNrOfBasisFunctions();++j){
	    element->basisFunctionCurl(j,point,curl);
	    curls.push_back(curl);
	}
	curlCache_[point]=curls;
    }
}

std::map<PointElementReferenceT,std::vector<NumericalVector> > FunctionCache::valueCache_;
std::map<PointElementReferenceT,std::vector<NumericalVector> > FunctionCache::curlCache_;

ElementInfos::ElementInfos(const Base::Element& element):inverse_(3,3),Jacobian_(3,3){
    PointElementReferenceT p(3);
    Geometry::PointPhysical pPhys(3);
    element.getReferenceGeometry()->getCenter(p);
    element.referenceToPhysical(p,pPhys);
    //not quite sure about the best way to implement this; this works for now
    if((pPhys[1]-0.5)*(pPhys[1]-0.5)+(pPhys[2]-0.5)*(pPhys[2]-0.5)<.25*.25){
    ///\bug Does not check that there are element boundaries at any of the discontinuities
    //if((pPhys[0]<0.3)||pPhys[0]>0.7||pPhys[1]<0.3||pPhys[1]>0.7){
	epsilon_=1;
    }else{
	epsilon_=1;
    }
    //the jacobian of a tetrahedron is constant.
    element.calcJacobian(p,Jacobian_);
    determinant_=Jacobian_.determinant();
    InvertAndTranspose(Jacobian_,inverse_);
    //cout<<"element "<<element.getID()<<" has jacobean determinant "<<determinant_<<endl;
}

void ElementInfos::makeFunctionValuesVector(const Base::Element* element, const PointElementReferenceT& point, std::vector< NumericalVector >& values){
    FunctionCache::getFunctionValuesVector(element,point,values);
    for(int j=0;j<element->getNrOfBasisFunctions();++j){
	//3D coordinate transformations -- something for a GPU?
	values[j]=inverse_*values[j];
    }        
}

void ElementInfos::makeFunctionCurlsVector(const Base::Element* element, const PointElementReferenceT& point, std::vector< NumericalVector >& curls){
    FunctionCache::getFunctionCurlsVector(element,point,curls);
    for(int j=0;j<element->getNrOfBasisFunctions();++j){
	//3D coordinate transformations -- something for a GPU?
	curls[j]=Jacobian_*curls[j]/determinant_;
    } 
}

MaxwellData::MaxwellData(int numberOfIntervals, int polynomialOrder):Sigma_(0),
    StabCoeff_(1.4),
    StartTime_(0),
    EndTime_(0.1),
    NumberOfIntervals_(numberOfIntervals),
    PolynomialOrder_(polynomialOrder){}
