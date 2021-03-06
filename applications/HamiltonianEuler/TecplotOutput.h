//
//  TecplotOutput.h
//
//
//  Created by Shavarsh Nurijanyan on 7/4/13.
//
//

#ifndef HPGEM_APP_TECPLOTOUTPUT_H
#define HPGEM_APP_TECPLOTOUTPUT_H
// To declare Physical and reference space points
#include "Geometry/PointPhysical.h"
#include "Geometry/PointReference.h"
#include "LinearAlgebra/MiddleSizeVector.h"
#include "Output/TecplotSingleElementWriter.h"
#include "InitialConditions.h"

using Geometry::PointPhysical;
using Geometry::PointReference;
#include <iomanip>
#include <iostream>
#include <fstream>

using namespace hpgem;

using std::ifstream;
using std::istream;
using std::ofstream;
using std::ostream;

class TecplotWriteFunction : public Output::TecplotSingleElementWriter<3> {

    // 	double exactSolutionP(const PhysSpacePoint<dim>& pPhys)
    // 	{
    // //std::cout<<"currentTime"<<CurrentTime<<endl;
    // 		return 2*std::cos(2*pi*(pPhys[0] + pPhys[1] + pPhys[2]) +
    // std::sqrt(3)*currTime/3)/(4*pi*pi);
    // 	}
    //
    // 	double exactSolutionU(const PhysSpacePoint<dim>& pPhys)
    // 	{
    // 		return (std::sqrt(3)*std::cos(2*pi*(pPhys[0] + pPhys[1] +
    // pPhys[2])
    // + std::sqrt(3)*currTime/3)+3*sin(2*pi*(pPhys[0] + pPhys[1] + pPhys[2])  +
    // std::sqrt(3)*currTime/3))/pi/2;
    // 	}
    //
    // 	double exactSolutionV(const PhysSpacePoint<dim>& pPhys)
    // 	{
    // 		return (std::sqrt(3)*std::cos(2*pi*(pPhys[0] + pPhys[1] +
    // pPhys[2])
    // + std::sqrt(3)*currTime/3)-3*sin(2*pi*(pPhys[0] + pPhys[1] + pPhys[2])  +
    // std::sqrt(3)*currTime/3))/pi/2;
    // 	}
    //
    // 	double exactSolutionW(const PhysSpacePoint<dim>& pPhys)
    // 	{
    // 		return -(2*std::sqrt(3)*std::cos(2*pi*(pPhys[0] + pPhys[1] +
    // pPhys[2])
    // + std::sqrt(3)*currTime/3))/pi/2;
    // 	}

   public:
    TecplotWriteFunction(ofstream* energyFile, ofstream* divFile,
                         ofstream* energyExfile, const ExactSolutionBase* vel,
                         ofstream* file1Dx, ofstream* file1D,
                         ofstream* file1DEx, ofstream* l2File)
        : velocity_(vel),
          maxError_(0),
          maxErrorU_(0),
          file1Dx_(file1Dx),
          file1D_(file1D),
          file1DEx_(file1DEx),
          energyFile_(energyFile),
          divFile_(divFile),
          energyExfile_(energyExfile),
          l2errorfile_(l2File),
          maxErrorPoint_() {
        *energyFile_ << "Time"
                     << " , "
                     << "Energy" << endl;
        *energyExfile_ << "Time"
                       << " , "
                       << "EnergyEx" << endl;
        *divFile_ << "Time"
                  << " , "
                  << "DIV" << endl;
    }
    void writeToTecplotFile(const Base::Element* element,
                            const PointReference<3>& pRef, ostream& os) {
        double lambda, u, uExact, lambdaExact, v, w, energy, vExact, wExact;
        double energyExact, uError, vError, wError;
        LinearAlgebra::MiddleSizeVector sol = element->getSolution(0, pRef);
        //        data.extractState(p, s);
        //
        lambda = sol[3];
        u = sol[0];
        v = sol[1];
        w = sol[2];
        //        uError	    = data.getUError();
        //        vError	    = data.getVError();
        //        wError	    = data.getWError();
        energy = 0;
        //        energyExact = 0;//data.energyExact_;
        //
        PointPhysical<3> pPhys = element->referenceToPhysical(pRef);
        //
        uExact = velocity_->getU(pPhys, time_);  // exactSolutionU(pPhys);
        vExact = velocity_->getV(pPhys, time_);  // exactSolutionV(pPhys);
        wExact = velocity_->getW(pPhys, time_);  // exactSolutionW(pPhys);/
        lambdaExact =
            velocity_->getP(pPhys, time_);  // exactSolutionP(pPhys);//
        //
        //
        os << u << "\t" << uExact << "\t" << (u - uExact) << "\t" << v << "\t"
           << vExact << "\t" << (v - vExact) << "\t" << w << "\t" << wExact
           << "\t" << (w - wExact) << "\t" << lambda << "\t" << lambdaExact
           << "\t" << (lambda - lambdaExact) << "\t" << energy;
        //        os <<
        //        u<<"\t"<<uExact<<"\t"<<uError<<"\t"<<v<<"\t"<<vExact<<"\t"<<vError<<"\t"<<w<<"\t"<<wExact<<"\t"<<wError<<"\t"<<
        //        lambda << "\t"<< lambdaExact << "\t" <<
        //        (lambda-lambdaExact)<<"\t"<<energy; if (maxError_<
        //        std::abs(lambda-lambdaExact))
        //        {
        //            maxError_ = std::abs(lambda-lambdaExact);
        //            maxErrorPoint_= pPhys;
        //        }
        //
        //            //			if (data.isInternal())
        //        {
        //            if (maxErrorU_<std::abs(u-uExact))
        //                maxErrorU_ = std::abs(u-uExact);
        //            if (maxErrorV_<std::abs(v-vExact))
        //                maxErrorV_ = std::abs(v-vExact);
        //            if (maxErrorW_<std::abs(w-wExact))
        //                maxErrorW_ = std::abs(w-wExact);
        //        }
    }
    void writeDiv(double max) {
        *divFile_ << time_ << " " << max << endl;
        cout << "Max error in discrete divergence=" << max << endl;
    }
    void outputErrors() {
        cout << "Maximum error of Pressure is =" << maxError_ << " at point ("
             << maxErrorPoint_[0] << ", " << maxErrorPoint_[1] << ", "
             << maxErrorPoint_[2] << endl;
        cout << "Maximum error of U is =" << maxErrorU_ << endl;
        cout << "Maximum error of V is =" << maxErrorV_ << endl;
        cout << "Maximum error of W is =" << maxErrorW_ << endl;
        // maxError_=maxErrorU_=maxErrorV_=maxErrorW_=0;
    }

    void writeEnergy() {
        //        double totalEnergy=0;
        //        double totalEnergyEx=0;
        //        double totalPressure=0;
        //        double totalErrorU=0;
        //        double totalErrorV=0;
        //        double totalErrorW=0;
        //        double totalErrorP=0;
        //        for (ElItType it = mesh.pubElCo.begin(); it !=
        //        mesh.pubElCo.end(); ++it)
        //        {
        //            ElementData& elData = data[(*it)->id()];
        //
        //            totalEnergy 	+= 	elData.getEnergy();
        //            totalEnergyEx 	+= 	elData.getExactEnergy();
        //            totalPressure 	+= 	elData.getPressureIntegral();
        //
        //                //				if (elData.internal_)
        //            {
        //                totalErrorU		+=	elData.getUError();
        //                totalErrorV		+=	elData.getVError();
        //                totalErrorW		+=	elData.getWError();
        //                totalErrorP     += 	elData.getPError();
        //            }
        //        }
        //
        //        *energyFile_ << std::setprecision (6)<<time_ <<"
        //        "<<totalEnergy <<endl; *energyExfile_<< std::setprecision (6)
        //        <<time_ <<" "<< totalPressure <<endl; *l2errorfile_<<
        //        std::setprecision (4)<<time_<<"\t\t"<< std::setprecision
        //        (10)<< std::sqrt(totalErrorU) <<"\t\t"<<
        //        std::sqrt(totalErrorV)<<"\t\t"<<
        //        std::sqrt(totalErrorW)<<"\t\t"<< std::sqrt(totalErrorP)<<endl;
    }

   public:
    double time_;
    double maxError_;
    double maxErrorU_;
    double maxErrorV_;
    double maxErrorW_;

    PointPhysical<3> maxErrorPoint_;

    void setOutputFiles(ofstream* file1Dx, ofstream* file1D,
                        ofstream* file1DEx) {
        file1Dx_ = file1Dx;
        file1D_ = file1D;
        file1DEx_ = file1DEx;
    }

   private:
    const ExactSolutionBase* velocity_;
    ofstream* file1Dx_;
    ofstream* file1D_;
    ofstream* file1DEx_;

    ofstream* energyFile_;
    ofstream* divFile_;
    ofstream* energyExfile_;
    ofstream* l2errorfile_;
};

#endif  // HPGEM_APP_TECPLOTOUTPUT_H
