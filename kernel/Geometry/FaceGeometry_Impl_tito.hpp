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


#include "FaceGeometry.hpp"   //---------------------------- this must be commented out before compiling the codes!!!!!

namespace Geometry {

    template<unsigned int DIM>
    class FaceGeometry;

    template<unsigned int DIM>
    FaceGeometry<DIM>::FaceGeometry(ElementGeometryT*       ptrElemL,
                                    const LocalFaceNrType&  localFaceNumL,
                                    ElementGeometryT*       ptrElemR,
                                    const LocalFaceNrType&  localFaceNumR) :
        leftElementGeom_(ptrElemL),
        localFaceNumberLeft_(localFaceNumL),
        rightElementGeom_(ptrElemR),
        localFaceNumberRight_(localFaceNumR),
        faceType_(INTERNAL)
    {
            // ~OC~
            // All that remains to be done for the Face-class's data is to
            // find the mapping from one side of the face to the other.  The
            // body of the constructor is entirely for that:
        
            // This code was taken from DGFaceInternal.  declare and fill
            // two arrays with the global node numbers of the face in the
            // (possibly) different permutations in which they occur seen
            // from the two elements
        
//        typename PhysicalGeometryT::VectorOfPointIndexesT globalNodeNrsL;
//        typename PhysicalGeometryT::VectorOfPointIndexesT globalNodeNrsR;
//            //
        
       // cout << "In Face constructor=" << endl;
//        cout << "Element left=" << *ptrElemL <<endl;
//        cout << "Element right=" << *ptrElemR << endl;
//        cout << "Local left face number=" <<localFaceNumL<<endl;
//        cout << "Local right face number=" <<localFaceNumR<<endl;
//        
//        cout << "Element Left=" << ptrElemL->getReferenceGeometry() <<endl;
//        cout << "Element Left Co DIMension mapping=" << ptrElemL->getReferenceGeometry()->getCoDIM1ReferenceGeometry(localFaceNumberLeft_) <<endl;
 
        
        //leftElementGeom_->getReferenceGeometry()->getCoDIM1ReferenceGeometry(localFaceNumberLeft_);
        
        std::vector<unsigned int> globalNodeNrsL;
        std::vector<unsigned int> globalNodeNrsR;
        
            // We try element left first because its the one we are sure exists.
        leftElementGeom_->getPhysicalGeometry()->getGlobalFaceNodeIndices(localFaceNumberLeft_, globalNodeNrsL);
        
            // This means internal face in old hpGEM (only left element for non-boundary)

        if (rightElementGeom_ != NULL && faceType_==INTERNAL)
        {
            rightElementGeom_->getPhysicalGeometry()->
            getGlobalFaceNodeIndices(localFaceNumberRight_,globalNodeNrsR);
        }
        else
        {
                // ~OC~ Otherwise just duplicate, the trafo will not be used
                // \bug: WTF is trafo?
                // (S.N) I bet, trafo is something connected to peridocity. This line will cause a creation of a boundary which has equal sets of PhysicalPoints (permutation included) from left and right. Thus the coDIMmappings are set accordingly.
                //
            globalNodeNrsR = globalNodeNrsL;
        }
            //  ~OC~
            //  Now we have a look at the global node indices.  The only use
            //  that is made of physical space information by the face is to
            //  generate the RefFace2RefFaceMapping.  If the two
            //  FaceDescriptors gave different sets of these indices then
            //  this is becoming a periodic face on one of the
            //  boundaries. In case the
            //  sets are identical, the ReferenceGeometry determines from
            //  the permutation of the index vectors which
            //  RefFace2RefFaceMapping needs to be used.
        
            //: The next to lines creating 2 std::sets, filling them with globalNodesL and globalNodesR. The sets are created via ctr which takes a to pointers and insert the data inbetween.
        cout<<"Left {";
        for(unsigned int i=0; i<globalNodeNrsL.size();++i)
        {
            cout<< globalNodeNrsL[i]<<" ";
        }
        cout<<"Left }"<<endl;
        
        cout<<"Right {";
        for(unsigned int i=0; i<globalNodeNrsR.size();++i)
        {
            cout<< globalNodeNrsR[i]<<" ";
        }
        cout<<"Right }"<<endl;
        SetOfGlobalNodes sL(globalNodeNrsL.data(), globalNodeNrsL.data() + globalNodeNrsL.size());
        SetOfGlobalNodes sR(globalNodeNrsR.data(), globalNodeNrsR.data() + globalNodeNrsR.size());

        
        if (sL == sR)// left set of global nodes for the face identical to right one. No periodicity detected!
        {
                // ~OC~
                // Equal sets algorithm checks all permutation too;
                // Find out which face2face-mapping index is needed
         //   cout << "The reference geometry " << endl;
         //   cout << this->referenceGeometryPtr() << endl;
            faceToFaceMapIndex_ = this->referenceGeometryPtr()->
            getCoDIM0MappingIndex(globalNodeNrsL, globalNodeNrsR);
        }
        else
        {
//-MTJ-start--------------
#ifdef MTJ
            //******************************
            // either the connected elements are refined or this is a periodic BC face
            //******************************

            // get intersection of left and right nodes set
            SetOfGlobalNodes commNodes;
            std::set_intersection( sL.begin(), sL.end(), sR.begin(), sR.end(), std::inserter( commNodes, commNodes.begin() ) );

            if (commNodes.size() > 0)
            // the connected elements share some nodes
            {
                SetOfGlobalNodes diffNodesL;
                std::set_difference( sL.begin(), sL.end(), sR.begin(), sR.end(), std::inserter(diffNodesL, diffNodesL.end() ) );

                SetOfGlobalNodes diffNodesR;
                std::set_difference( sR.begin(), sR.end(), sL.begin(), sL.end(), std::inserter(diffNodesR, diffNodesR.end() ) );

    //             typedef typename Geometry::PhysicalGeometry<DIM>::NodeContainerPtrType NodeContainerPtrType ;
    //             NodeContainerPtrType PCPtr = _elL->physicalGeometry()->getPointContainerPtr();

                bool foundCollinear = false;
                for(unsigned int i=0; i<globalNodeNrsL.size(); ++i)
                {
                  PointPhysical<DIM> ppR;
    //               PCPtr->getPoint(globalNodeNrsR(i), ppR);
                  ppR = *(rightElementGeom_->getPhysicalGeometry()->getNodePtr(globalNodeNrsR(i)));
                  
                  SetOfGlobalNodes::iterator itFound = commNodes.find(globalNodeNrsR(i));
                  if (itFound == commNodes.end())
                  {
                    // do collinear test for all possible pairs
                    for (SetOfGlobalNodes::iterator itCommon=commNodes.begin(); (itCommon!=commNodes.end()); ++itCommon)
                    {
                      PointPhysical<DIM> pp0;
    //                   PCPtr->getPoint(*itCommon, pp0);
                      pp0 = *(leftElementGeom_->getPhysicalGeometry()->getNodePtr(*itCommon));
                      
                      for (SetOfGlobalNodes::iterator itDiffL=diffNodesL.begin(); (itDiffL!=diffNodesL.end()); ++itDiffL)
                      {
                        if (diffNodesR.find(globalNodeNrsR(i)) == diffNodesR.end()) continue;
                        
                        PointPhysical<DIM> ppL;
    //                     PCPtr->getPoint(*itDiffL, ppL);
                        ppL = *(leftElementGeom_->getPhysicalGeometry()->getNodePtr(*itDiffL));
                        
                        PointPhysical<DIM> dL(ppL-pp0);
                        PointPhysical<DIM> dR(ppR-pp0);

                        double ratio = 0.;
                        unsigned int d;
                        for (d=0; d<DIM; ++d)
                        {
                          if (std::abs(dL[d])>SmallerDoubleThanMinimalSizeOfTheMesh)
                          {
                            ratio = dR[d]/dL[d];
                            break;
                          }
                        }

                        bool isCollinear = true;
                        for (;d<DIM; ++d)
                        {
                          if (std::abs(dL[d])>SmallerDoubleThanMinimalSizeOfTheMesh)
                          {
                            if (std::abs(dR[d]/dL[d] - ratio) > SmallerDoubleThanMinimalSizeOfTheMesh)
                            {
                              isCollinear = false;
                              break;
                            }
                          }
                        }
                        
                        if (isCollinear) 
                        {
                          // replace nodes on the right list
                          diffNodesR.erase(diffNodesR.find(globalNodeNrsR(i)));
                          globalNodeNrsR(i) = *itDiffL;
                          foundCollinear = true;
                        }
                      }  // for itDiffL
                    }  // for itCommon
                  }  // if itFound
                }  // for i
                
                if (!foundCollinear) throw "Face: incorrect pairs!";
            } // end if commNodes.size() > 0

            else
            {
#endif
//-MTJ-end--------------
              
                    // ~OC~
                VectorOfLocalNodes localNodeNrsL;
                VectorOfLocalNodes localNodeNrsR;
                leftElementGeom_->getPhysicalGeometry()->getLocalFaceNodeIndices(localFaceNumberLeft_, localNodeNrsL);
                rightElementGeom_->getPhysicalGeometry()->getLocalFaceNodeIndices(localFaceNumberRight_, localNodeNrsR);
                
                
                unsigned int periodicDim = Geometry::MaxInteger; // initialize to large value
                
                PointPhysicalT pointPhysical;
                
                std::vector<PointPhysicalT> projectedPointsL; // create empty!
                std::vector<PointPhysicalT> projectedPointsR;
                
                for (unsigned int i = 0; i < localNodeNrsL.size(); ++i)
                {
                    leftElementGeom_->getPhysicalGeometry()->getNodeCoordinates(localNodeNrsL[i], pointPhysical);
                    
                    projectedPointsL.push_back(pointPhysical);
                    
                    rightElementGeom_->getPhysicalGeometry()->getNodeCoordinates(localNodeNrsR[i], pointPhysical);
                    
                    projectedPointsR.push_back(pointPhysical);
                }
                
                    // We want to find which direction is the periodic one.
                unsigned int test[DIM];
                
                for (unsigned int d = 0; d < DIM; ++d) { test[d] = 0; }
                
                pointPhysical = projectedPointsL[0]; //Find out which one is peridoci by comparing all coordinates of one point to other points, for all DIMension. If DIMension is the same for every Point, than it is a face that DIMension in question.
                for (unsigned int i = 0; i < localNodeNrsL.size(); ++i)
                {
                    for (unsigned int d = 0; d < DIM; ++d)
                    {
                        if (std::abs(projectedPointsL[i][d]-pointPhysical[d]) < SmallerDoubleThanMinimalSizeOfTheMesh)
                        {
                            test[d] += 1;
                        }
                    }
                }
                    //check if for that DIMension is fullhouse
                for (unsigned int d = 0; d < DIM; ++d)
                {
                    if (test[d] == localNodeNrsL.size())
                    {
                        periodicDim = d;
                    }
                }
                    //  periodicDim should be our DIMension, otherwise error!
                if (periodicDim==Geometry::MaxInteger)
                    cout << "Shit happened. Trained monkeys are on the way to fix the problem, probably." << std::endl;
                
                
                PhysicalPointOnTheFaceT ppL;
                PhysicalPointOnTheFaceT ppR;
                
                    // next we eliminate the periodic DIMension from a Point coordinates & a
                for(unsigned int i = 0; i < globalNodeNrsR.size(); ++i)
                {
                    for (unsigned int d = 0; d < periodicDim; ++d)
                    {
                        
                        projectedPointsR[i][d];
                    }
                    for (unsigned int d = periodicDim+1; d < DIM; ++d)
                    {
                        ppR[d-1] = projectedPointsR[i][d];
                    }
                    
                    for (unsigned int j = 0; j < globalNodeNrsL.size(); ++j)
                    {
                        for (unsigned int d = 0; d < periodicDim; ++d)
                        {
                            ppL[d] = projectedPointsL[j][d];
                        }
                        for (unsigned int d = periodicDim+1; d < DIM; ++d)
                        {
                            ppL[d-1] = projectedPointsL[j][d];
                        }
                        if (ppL == ppR)//check if now they are the same, then overwrite the Right nodeList, with LeftNodeList
                        {
                            globalNodeNrsR[i] = globalNodeNrsL[j];
                        }
                    }
                }
//-MTJ-start--------------
#ifdef MTJ
            }
#endif
//-MTJ-end--------------
            
            faceToFaceMapIndex_ = this->referenceGeometryPtr()->getCoDIM0MappingIndex(globalNodeNrsL, globalNodeNrsR);
        }
    }

        //! Ctor for boundary faces.
	template<unsigned int DIM>
    FaceGeometry<DIM>::FaceGeometry(ElementGeometryT* ptrElemL, const LocalFaceNrType& localFaceNumL, const FaceType& boundaryLabel):
    leftElementGeom_(ptrElemL),
    rightElementGeom_(NULL),
    localFaceNumberLeft_(localFaceNumL),
    localFaceNumberRight_(Geometry::MaxInteger),
    faceToFaceMapIndex_(0),
    faceType_(boundaryLabel)
    { }
        /// The referenceGeometry is returned. It is taken from left element, to always ensure its existence
        // should use the same interface as in the Element!
    template<unsigned int DIM>
    const ReferenceGeometry<DIM-1>*
    FaceGeometry<DIM>::referenceGeometryPtr()const
    {
        return leftElementGeom_->getReferenceGeometry()->getCoDIM1ReferenceGeometry(localFaceNumberLeft_);
    }
    
    template<unsigned int DIM>
    const ReferenceGeometry<DIM-1>*
    FaceGeometry<DIM>::getReferenceGeometry()const
    {
        return leftElementGeom_->getReferenceGeometry()->getCoDIM1ReferenceGeometry(localFaceNumberLeft_);
    }
 
	
//        //! Return the unique number of this face.
//	const IDType& id() const
//    {
//		return _id;
//    }
    
 

		
	/*! \brief Return a wrapped elemental function which can be evaluated at
	 *  reference face points. In this case the function is defined on the
	 *  left element. */
//	template <class FType>
//	ElementDataTraceL<DIM, FType, Face>
//	transform2RefFaceL(FType f) const
//    {
//		return ElementDataTraceL<DIM, FType, Face>(//this,
//                                                   f);
//    }
	
	/*! \brief Return a wrapped elemental function which can be evaluated at
	 *  reference face points. In this case the function is defined on the
	 *  right element. */
//	template <class FType>
//	ElementDataTraceR<DIM, FType, Face>
//	transform2RefFaceR(FType f) const
//    {
//		return ElementDataTraceR<DIM, FType, Face>(//this,
//                                                   f);
//    }
	
	/*! Map a point in coordinates of the reference geometry of the face to
	 *  the reference geometry of the left (L) element. */
    template<unsigned int DIM>
    void
    FaceGeometry<DIM>::mapRefFaceToRefElemL(const ReferencePointOnTheFaceT& pRefFace,
                                     ReferencePointT& pRefEl) const
    {
		leftElementGeom_->getReferenceGeometry()->getCoDIM1MappingPtr(localFaceNumberLeft_)->transform(pRefFace, pRefEl);
    }
	
	/*! Map a point in coordinates of the reference geometry of the face to
	 *  the reference geometry of the right (R) element. */
    template<unsigned int DIM>
    void
    FaceGeometry<DIM>::mapRefFaceToRefElemR(const ReferencePointOnTheFaceT& pRefFace,
                                     ReferencePointT& pRefEl) const
    {
            // In the L function we have assumed the point pRefFace to be
            // given in coordinates of the system used by the reference face
            // on the L side. That means that now to make sure that the
            // point transformed from the right side of the face onto the
            // right element is the same as the one on the left side of the
            // face, we have to use the refFace2RefFace mapping.
		ReferencePointOnTheFaceT pOtherSide;
 		mapRefFaceToRefFace(pRefFace, pOtherSide);
		rightElementGeom_->getReferenceGeometry()->getCoDIM1MappingPtr(localFaceNumberRight_)->transform(pOtherSide, pRefEl);
    }
	
	/*! Map from reference face coordinates on the left side to those on the
	 *  right side. */
	template<unsigned int DIM>
    void
    FaceGeometry<DIM>::mapRefFaceToRefFace(const ReferencePointOnTheFaceT& pIn,
                                     ReferencePointOnTheFaceT& pOut) const
    {
		referenceGeometryPtr()->getCoDIM0MappingPtr(faceToFaceMapIndex_)->transform(pIn, pOut);
    }
	
//	typedef auto_ptr<const Ref2RefSpaceMapping<DIM-1, DIM> >
//	RefFace2RefElementMapping;
	
        //! Return a Mapping (not pointer or reference! Ok, wrapped by auto_ptr)
	template<unsigned int DIM>
    typename FaceGeometry<DIM>::RefFaceToRefElementMapping
    FaceGeometry<DIM>::refFaceToRefElemMapL() const
    {
		return
        RefFaceToRefElementMapping(leftElementGeom_->getReferenceGeometry()->getCoDIM1MappingPtr(localFaceNumberLeft_));
    }
	
        //! Return a mapping to the right reference element.
	template<unsigned int DIM>
    typename FaceGeometry<DIM>::RefFaceToRefElementMapping 
    FaceGeometry<DIM>::refFaceToRefElemMapR() const
    {
		const MappingReferenceToReference<DIM-1, DIM-1>* const m1Ptr = this->referenceGeometryPtr()->getCoDIM0MappingPtr(faceToFaceMapIndex_);
		const MappingReferenceToReference<DIM-1, DIM>* const m2Ptr =  leftElementGeom_->getReferenceGeometry()->getCoDIM1MappingPtr(localFaceNumberRight_);
		
		return RefFaceToRefElementMapping(new ConcatenatedMapping<DIM-1, DIM-1, DIM>(*m1Ptr, *m2Ptr));
    }
    
    /*!Compute the normal vector at a given point (in face coords).
     
     Faces turn up in DG discretizations as domains for some of the
     integrals. For these integrals, the face must be able to provide the
     geometric information. This is completely included in the normal
     vector, which gives
     <UL>
     <LI> the direction of the normal vector oriented from left (L) to
     right (R) element; that way, for a boundary face, it is an external
     normal vector;
     <LI> the transformation of the integration element between reference
     face and physical space: this scalar is given by the 2-norm of the
     returned normal vector.
     </UL> */
    template<unsigned int DIM>
    void
    FaceGeometry<DIM>::getNormalVector(const ReferencePointOnTheFaceT& pRefFace, PointPhysicalT& v) const
    {
        	// first Jacobian (mapping reference face -> reference element)
		Jacobian<DIM-1, DIM>    j1;
		
		leftElementGeom_->getReferenceGeometry()->getCoDIM1MappingPtr(localFaceNumberLeft_) // this is the refFace2RefElemMap
        ->calcJacobian(pRefFace, j1);
        
            // second Jacobian (mapping reference element -> phys. element),
            // for this we first need the points coordinates on the
            // reference element
		ReferencePointT pRefElement;
		mapRefFaceToRefElemL(pRefFace, pRefElement);
        
        Jacobian<DIM, DIM> j2;
        
		leftElementGeom_->calcJacobian(pRefElement, j2);
        
		Jacobian<DIM-1, DIM> j3;
        
        j2.multiplyJacobiansInto(j1, j3);
        
        j3.computeWedgeStuffVector(v);
        
        double det = j2.determinant();
		
        //sgn==(x > 0) - (x < 0)
        v *= ((det>0)-(det<0)) *
        OutwardNormalVectorSign<DIM>(leftElementGeom_->getReferenceGeometry()->getCoDIM1MappingPtr(localFaceNumberLeft_));
    }
};