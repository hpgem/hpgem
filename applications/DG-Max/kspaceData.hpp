#include <vector>
#include "LinearAlgebra/NumericalVector.hpp"

/**
 * The tesselation in k-space does not fit the hpGEM-framework because of its continuous nature
 * this is a specialized class to store the relevant information and provide the relevant computations
 */
class KspaceData{
  
    std::vector<LinearAlgebra::NumericalVector> kpoints_,deltak_;
    std::vector<std::vector<double> > omegaAtKpoints_;
    std::vector<std::vector<int> > elements_;
    std::vector<std::vector<LinearAlgebra::NumericalVector> > functionValuesAtKpoints_;
    int current_;
    int minimumsize_;
    
public:
  
    /**
     * constructor takes a desired number of k-point per direction and constructs a tesselation of small tetrahedra
     * that can be used for linear interpolation. the actual number of k-points where the eigenvalue problem needs to be solved
     * is (n)(n+1)(n+2)/6. makes sure the first k-point is the origin so no matrix updates are needed for the first pass
     */
    KspaceData(int pointsPerDirection);
    
    /**
     * returns true as long as there is still a point that needs computing
     */
    bool hasNextPoint();
    
    /**
     * iterator over the k-points. Call this when you are done for the current k-point and it gets you a update-vector
     */
    LinearAlgebra::NumericalVector& nextPoint();
    
    /**
     * sets the found values of omega at the current k-point. Make sure you set the same number of omega each k-point and that they are sorted in the same way (accounting for band crossings)
     */
    void setOmega(std::vector<double>& omega);
    
    /**
     * optional: sets the value of f from the expression int(f*delta(omega-omega_n))dk
     */
    void setFunctionValues(std::vector<LinearAlgebra::NumericalVector>& functionValues);
    
    /**
     * returns the density of states for a given value omega
     */
    void getIntegral(double omega, LinearAlgebra::NumericalVector& result);
};