//bool Base::MeshManipulator:readMesh(FILE file)
//{
    // TODO Implement reading mesh files.
    // For now, creating basic mesh for testing purposes
//    createElement(new ElementT())
//}


//This is a small struct used for storing the halfFaces that are used 
struct halfFaceDescription {
    unsigned int firstNode;
    unsigned int secondNode;
    unsigned int elementNum;
    unsigned int localFaceIndex;
};
 ///Half faces are used by the faceFactory. These are face which only know one element they are associated with
typedef struct halfFaceDescription halfFaceDescriptionT;

/// \brief This is the function that checks if two halfFaces nned to swapped or not.
/// \details 
/*!
 *  First checks the firstNode index and if this is the same, swap based on the secondNode index.
 *  Therefore in the end that will be sorted based on firstNode and secondly on secondNode.
 !*/
bool compareHalfFace(halfFaceDescriptionT first, halfFaceDescriptionT second)
{
    if (first.firstNode<second.firstNode) return true;
    
    if ((first.secondNode<second.secondNode) && (first.firstNode==second.firstNode)) return true;
    
    return false;
}

template<>
void
MeshManipulator<1>::createBasisFunctions(unsigned int order)
{
    Base::BasisFunctionSet<1>* bFset1 = new Base::BasisFunctionSet<1>(order);
    
    Base::AssembleBasisFunctionSet_1D_Ord3_A0(*bFset1);
    
    collBasisFSet_.push_back(bFset1);
}

template<>
void
MeshManipulator<2>::createBasisFunctions(unsigned int order)
{
    Base::BasisFunctionSet<2>* bFset1 = new Base::BasisFunctionSet<2>(order);
    
    Base::AssembleBasisFunctionSet_2D_Ord2_A1(*bFset1);
    
    collBasisFSet_.push_back(bFset1);
}

template<>
void
MeshManipulator<3>::createBasisFunctions(unsigned int order)
{
    Base::BasisFunctionSet<3>* bFset1 = new Base::BasisFunctionSet<3>(order);
    
    Base::AssembleBasisFunctionSet_3D_Ord2_A1(*bFset1);
    
    collBasisFSet_.push_back(bFset1);
}

template<unsigned int DIM>
void MeshManipulator<DIM>::createRectangularMesh(PointPhysicalT BottomLeft, PointPhysicalT TopRight, const  std::vector<unsigned int> linearNoElements)
{
    if (linearNoElements.size() != DIM)
    {
        cout << "The number of Linear Intervals has to map the size of the problem and current it does not"<<endl;
        throw(10);
    }
    
    //Stage 1 : Precompute some required values;
    ///////
    
    //This store the size length of the domain i.e. it is DIM sized vector
    PointPhysicalT delta_x(DIM);
    
    for (int i=0;i<DIM;i++)
    {delta_x[i]=(TopRight[i]-BottomLeft[i])/(linearNoElements[i]);}

    

    
    //This stores the number of nodes in each coDIMension i.e. if you have 2 by 2 element it is 3 nodes 
    std::vector<int> numOfNodesInEachSubspace(DIM), numOfElementsInEachSubspace(DIM);
    
    numOfNodesInEachSubspace[0]=1;
    numOfElementsInEachSubspace[0]=1;
    
    //This will be the total number of nodes required in the problem
    unsigned int totalNumOfNodes,totalNumOfElements, verticesPerElement;
    
    totalNumOfNodes=(linearNoElements[0]+1);
    totalNumOfElements=(linearNoElements[0]);
    verticesPerElement=2;
    int powerOf2;
    
    for (int iDIM=1; iDIM<DIM;++iDIM)
    {
        totalNumOfNodes*=(linearNoElements[iDIM]+1);
        totalNumOfElements*=(linearNoElements[iDIM]);
        verticesPerElement*=2;
        
        numOfElementsInEachSubspace[iDIM]=numOfElementsInEachSubspace[iDIM-1]*(linearNoElements[iDIM-1]);
        numOfNodesInEachSubspace[iDIM]=numOfNodesInEachSubspace[iDIM-1]*(linearNoElements[iDIM-1]+1);
        
    }
    
    //temp point for storing the node locations
    PointPhysicalT x;
    
    
    ///Stage 2 : Create the nodes
    //Now loop over all the nodes and calculate the coordinates for rach DIMension (this makes the algorithm independent of DIMension
    for (int nodeIndex=0; nodeIndex<totalNumOfNodes; ++nodeIndex)
    {
        int nodeIndexRemain=nodeIndex;
        
        
        
        for (int iDIM=DIM-1;iDIM>-1;--iDIM)
        {
            x[iDIM] = BottomLeft[iDIM] + (nodeIndexRemain/numOfNodesInEachSubspace[iDIM]*delta_x[iDIM]);
            nodeIndexRemain %=numOfNodesInEachSubspace[iDIM];
        }
    
                                          
        //actally add the point
        points_.push_back(x);

    }
        
        
        
    //Stage 3 : Create the elements
   
    //Create an Instance of a referenceSquare
   
    
    //Geometry::ReferenceGeometry<DIM>* referenceGeometry;
    
    void* referenceGeometryVoidPtr;
    
    Geometry::ReferenceGeometry<DIM>* referenceGeometry;
    
    if (DIM==1)
    {
        
        referenceGeometryVoidPtr = &Geometry::ReferenceLine::Instance();
    }
    else if  (DIM==2)
    {
        referenceGeometryVoidPtr = &Geometry::ReferenceSquare::Instance();
    }
    else if (DIM==3)
    {
        referenceGeometryVoidPtr = &Geometry::ReferenceCube::Instance();
    }
    
    referenceGeometry=static_cast<Geometry::ReferenceGeometry<DIM>* > (referenceGeometryVoidPtr);
    
    std::vector<Base::Element<DIM>* > tempElementVector(totalNumOfElements);
    
    std::vector<unsigned int> elementNdId(DIM), vertexNdId(DIM), globalVertexID(verticesPerElement);
    
    //This will store a two array, which goes elementNumber and the n-DIMensional coordinate of that element. This is only used for the face creation (step 4).
    //std::vector<std::vector<unsigned int> > allElementNdId;
    //allElementNdId.resize(totalNumOfElements);
    
  
    //elementNdId is DIM coordinate of the bottom left node i.e. in two (0,0), (1,0) ,(2,0) ... etc are the first three (if at least three elements in x)
    for (int elementIndex=0; elementIndex<totalNumOfElements; ++elementIndex)
    {
        int numElementsRemaining=elementIndex;
        
        for (int iDIM=DIM-1; iDIM>-1;--iDIM)
        {
            elementNdId[iDIM]=numElementsRemaining/numOfElementsInEachSubspace[iDIM];
            numElementsRemaining %= numOfElementsInEachSubspace[iDIM];
        }
    //    allElementNdId[elementIndex].resize(DIM);
    //    allElementNdId[elementIndex]=elementNdId;
        
        // vertexNdId are the DIM coordinate of each vertex in the element with vertexNdId[0] being the bottom left
        for (int i=0; i<verticesPerElement; ++i)
        {
            powerOf2 = 1;
            for (int iDIM=0; iDIM<DIM; ++iDIM)
            {
                vertexNdId[iDIM] = elementNdId[iDIM] + ((i & powerOf2) !=0);
                powerOf2 *=2;
            }
            globalVertexID[i]=vertexNdId[0];
                
            //Now map to the one DIMensional global ID
            for (int iDIM=1;iDIM<DIM;++iDIM){globalVertexID[i] += vertexNdId[iDIM]*numOfNodesInEachSubspace[iDIM];}
        }
            
        tempElementVector[elementIndex]=addElement(globalVertexID,referenceGeometry);
    }
    
    
    //Stage 4 : Create the faces
    /// \bug Only implemented for 1D and 2D at the momennt
    // Note the linear number of faces in each direction is linearNoElements-1;
    
    //loop over all the elements
    //for (int elementIndex=0;elementIndex<totalNumOfElements;elementIndex++)
//    {
//        //loop over the DIMensions of the element. Element create internal faces in possitive direction and boundary faces
//        for (int DIMIndex=0; DIMIndex<DIM;DIMIndex++)
//        {
//            //Left x boundary
//            if (allElementNdId[elementIndex][DIMIndex]==0)
//            {
//                addFace(tempElementVector[elementIndex],
//            }
//        
//        }
//    }
    
   // rectangularCreateFaces1D(&tempElementVector, &linearNoElements);
    switch (DIM)
    {
        case 1:
            rectangularCreateFaces1D(tempElementVector, linearNoElements);
            break;
            
        case 2:
            rectangularCreateFaces2D(tempElementVector, linearNoElements);
            break;
            
        case 3:
            rectangularCreateFaces3D(tempElementVector, linearNoElements);
            break;
       
        default:
            throw("Face generator not implemented in this DIMension");
            
    }//end case
    
   
}


template<unsigned int DIM>
void MeshManipulator<DIM>::rectangularCreateFaces1D(std::vector<Base::Element<DIM>* >& tempElementVector, const std::vector<unsigned int>& linearNoElements)
{
    
    //First the internal faces
    for (int i=0;i<linearNoElements[0]-1;i++)
    {
        addFace(tempElementVector[i],1,tempElementVector[i+1],0);
    }
    
    
    if (periodicX_==1)
        
    {
        addFace(tempElementVector[0],0,tempElementVector[linearNoElements[0]-1],1);
    }   
    else 
    {
        addFace(tempElementVector[0],0,NULL,0);
        addFace(tempElementVector[linearNoElements[0]-1],1,NULL,0);
    }
    
}


template<unsigned int DIM>
void MeshManipulator<DIM>::rectangularCreateFaces2D(std::vector<Base::Element<DIM>* >& tempElementVector, const std::vector<unsigned int>& linearNoElements)
{
    
    //first do the faces with normals pointing in the x-direction
    Geometry::FaceType xFace = Geometry::WALL_BC;
    unsigned int index;
    
    for (int j=0; j<linearNoElements[1]; j++)
    {
        for (int i=0; i<linearNoElements[0]-1;i++)
        {
            
            index=i+j*linearNoElements[0];
            addFace(tempElementVector[index],2,tempElementVector[index+1],1);
            
        }
        
        
        index=j*linearNoElements[0];
        if (periodicX_==1)
        {
            addFace(tempElementVector[index],1,tempElementVector[index+linearNoElements[0]-1],2);
        }
        else 
        {
            
            //Left bounday face
            
            cout<<"index1="<<index<<endl;
            addFace(tempElementVector[index],1,NULL,0);
            
            //Right boundary face
            index =index+linearNoElements[0]-1;
            cout<<"index2="<<index<<endl;
            addFace(tempElementVector[index],2,NULL,0);
        }
    }
    
    //now do the faces with normals pointing in the y-direction
    Geometry::FaceType yFace=Geometry::WALL_BC;
    for (int j=0; j<linearNoElements[0]; j++)
    {
        for (int i=0; i<linearNoElements[1]-1;i++)
        {
            index=j+i*linearNoElements[0];
            cout<<"index3="<<index<<endl;
            cout<<"index3.5="<<index+linearNoElements[0]<<endl;
            
            addFace(tempElementVector[index],0,tempElementVector[index+linearNoElements[0]],3);
            
        }
        
        index=j;
        
        if (periodicY_==1)
        {
            
            addFace(tempElementVector[index],0,tempElementVector[index+(linearNoElements[1]-1)*linearNoElements[0]],3);
            
        }
        else 
        {
            
            //Bottom boundary face
            
            cout<<"index4="<<index<<endl;
            addFace(tempElementVector[index],0,NULL,0);
            
            //Top boundary face
            index =index+(linearNoElements[1]-1)*linearNoElements[0];
            
            cout<<"index5="<<index<<endl;
            addFace(tempElementVector[index],3,NULL,0);
        }
        
    }
    
    
}

template<unsigned int DIM>
void MeshManipulator<DIM>::rectangularCreateFaces3D(std::vector<Base::Element<DIM>* >& tempElementVector, const std::vector<unsigned int>& linearNoElements)
{
    
    unsigned int index;
    //first do the faces in x-direction
    //counter in z
    for(int k=0;k<linearNoElements[2];k++)
    {
        //counter in y
        for (int j=0; j<linearNoElements[1]; j++)
        {
            //counter in x
            for (int i=0; i<linearNoElements[0]-1;i++)
            {
            
                index=i+j*linearNoElements[0]+k*linearNoElements[0]*linearNoElements[1];
                addFace(tempElementVector[index],3,tempElementVector[index+1],2);
            
            }//end loop over x
        
        
            index=j*linearNoElements[0]+k*linearNoElements[0]*linearNoElements[1];
            if (periodicX_==1)
            {
                addFace(tempElementVector[index],2,tempElementVector[index+linearNoElements[0]-1],3);
            }
            else 
            {
            
                //Left boundary
            
                cout<<"index1="<<index<<endl;
                addFace(tempElementVector[index],2,NULL,0);
            
                //Right boundary face
                index =index+linearNoElements[0]-1;
                cout<<"index2="<<index<<endl;
                addFace(tempElementVector[index],3,NULL,0);
            } // end boundary cases
            
        } // end loop over y
    } //end loop over z
    
    
    //Now do the faces pointing in the y-direction
    //count in z
    for(int k=0;k<linearNoElements[2];k++)
    {
        //counter in x
        for (int j=0; j<linearNoElements[0]; j++)
        {
            //counter in y
            for (int i=0; i<linearNoElements[1]-1;i++)
            {
                
                index=j+i*linearNoElements[0]+k*linearNoElements[0]*linearNoElements[1];
                addFace(tempElementVector[index],4,tempElementVector[index+linearNoElements[0]],1);
                
            }//end loop over x
            
            
            index=j+k*linearNoElements[0]*linearNoElements[1];
            if (periodicY_==1)
            {
                cout << "Hmmm this is the problem " << endl;
                cout << index << " , " << index+(linearNoElements[1]-1)*linearNoElements[0] << endl;
                addFace(tempElementVector[index],1,tempElementVector[(index+(linearNoElements[1]-1)*linearNoElements[0])],4);
            }
            else 
            {
                
                //front bounday face
                
                cout<<"index1="<<index<<endl;
                addFace(tempElementVector[index],1,NULL,0);
                
                //Back boundary face
                index =index+(linearNoElements[1]-1)*linearNoElements[0];
                cout<<"index2="<<index<<endl;
                addFace(tempElementVector[index],4,NULL,0);
            } // end boundary cases
            
        } // end loop over y
    } //end loop over z
    
        
    
    //Now finally do the face in the z-direction
    
    //count in x direction 
    for(int k=0;k<linearNoElements[0];k++)
    {
        //counter in y
        for (int j=0; j<linearNoElements[1]; j++)
        {
            //counter in z
            for (int i=0; i<linearNoElements[2]-1;i++)
            {
                
                index=k+j*linearNoElements[0]+i*linearNoElements[0]*linearNoElements[1];
                addFace(tempElementVector[index],5,tempElementVector[index+linearNoElements[0]*linearNoElements[1]],0);
                
            }//end loop over x
            
            
            index=k+j*linearNoElements[0];
            if (periodicZ_==1)
            {
                cout << "This is the final problem " <<endl;
                cout << index << " , ";
                cout << index+(linearNoElements[2]-1)*linearNoElements[1]*linearNoElements[0] << endl;
                addFace(tempElementVector[index],0,tempElementVector[index+(linearNoElements[2]-1)*linearNoElements[1]*linearNoElements[0]],5);
            }
            else 
            {
                
                //bottom boundary
                cout<<"index1="<<index<<endl;
                addFace(tempElementVector[index],0,NULL,0);
                
                //Top boundary face
                index =index+(linearNoElements[2]-1)*linearNoElements[1]*linearNoElements[0];
                cout<<"index2="<<index<<endl;
                addFace(tempElementVector[index],5,NULL,0);
            } // end boundary cases
            
        } // end loop over y
    } //end loop over z
    
}

template<unsigned int DIM>
void MeshManipulator<DIM>::readCentaurMesh(std::string filename)
{

    //First open the file
    std::ifstream centaurFile;
    
    centaurFile.open(filename.c_str(), std::ios::binary);
	if (!centaurFile.is_open())
	{
	    throw("Cannot open Centaur meshfile.");
	}
    
    switch (DIM)
    {
        case 2: 
            readCentaurMesh2D(centaurFile);
            break;
        default:
            std::cerr<<"Centaur mesh reader has not been implemented in this DIMension" << std::endl;
    }
    
    
    //Finally close the file
    centaurFile.close();
    
}

template<unsigned int DIM>
void MeshManipulator<DIM>::readCentaurMesh2D(std::ifstream& centaurFile)
{
  
    
    //These are used to check the length of the read lines to check for read errors
    int sizeOfLine;	   
    
    //This first value in the centaur file is the size of each line in the file;
    centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
    
    // Version number of the Centaur mesh file 0.1 if using Tito's matlab generator
    float version;		      
	centaurFile.read(reinterpret_cast<char*>(&version), sizeof(version));
    std::cout << "This read mesh is in Centaur version " << version << " format" <<std::endl;
    
    // Centaur File Type <0 is two DIMensional and >0 is three DIMensional
    int centaurFileType;			
	centaurFile.read(reinterpret_cast<char*>(&centaurFileType),sizeof(centaurFileType));
    
    
    
    if (centaurFileType<0)
    {
        std::cout << "Reading a two DIMensional centaur mesh" <<std::endl;
        
        //The rest of the first line is junk
        char junk[1024];
        
        int checkInt;
        centaurFile.read(&junk[0], sizeOfLine-sizeof(version)-sizeof(centaurFileType));
        
        //Check the first line was read correctly : each line in centaur start and end with the line size as a check
        centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
        if (checkInt!=sizeOfLine) {std::cerr << "Error in centaur file " << std::endl; return;}
        
        //Start the second line
        centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
        
        //Next read the total number of nodes
        int numberOfNodes;
	    centaurFile.read(reinterpret_cast<char*>(&numberOfNodes),sizeof(numberOfNodes));
        std::cout << "File contains " << numberOfNodes <<" nodes" <<std::endl;
        
        //Check the second line was read correctly : each line in centaur start and end with the line size as a check
        centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
        if (checkInt!=sizeOfLine) {std::cerr << "Error in centaur file " << std::endl; return;}
        
        
        //Now we will read in all the nodes
        centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
        double nodeCoord[2];
        PointPhysicalT nodeCoordPointFormat;
        for (int i=0; i<numberOfNodes; i++)
        {
            // Reads the x and y coordinates of each node.
            centaurFile.read(reinterpret_cast<char*>(&nodeCoord[0]), sizeof(nodeCoord));
            // pass the node to the nodelist.
            
            //Covert from *double to hpGEM PointPhysical format
            nodeCoordPointFormat[0]=nodeCoord[0];
            nodeCoordPointFormat[1]=nodeCoord[1];
            points_.push_back(nodeCoordPointFormat);
           
	    }
        //Now check the node line was read correctly : each line in centaur start and end with the line size as a check
        centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
        if (checkInt!=sizeOfLine) {std::cerr << "Error in centaur file " << std::endl; return;}
        
        
        
        //Now check how many triangle in the file
        centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
        // Number of triangular elements
        int numberOfTriangles;		
	    centaurFile.read(reinterpret_cast<char*>(&numberOfTriangles),sizeof(numberOfTriangles));
        std::cout << "File contains " <<numberOfTriangles<< " triangle(s)" << std::endl;
        
        //Check the line was read correctly : each line in centaur start and end with the line size as a check
        centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
        if (checkInt!=sizeOfLine) {std::cerr << "Error in centaur file " << std::endl; return;}
        
        centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
        /// \bug Need to add a triangle reader
        if (numberOfTriangles>0)
        {
            unsigned int triGlobalNodeIndexes[3];
            std::vector<unsigned int> globalNodeIndexes(3);
            Geometry::ReferenceGeometry<2>* referenceTriangle = &Geometry::ReferenceTriangle::Instance();
         
            for (int i=0; i<numberOfTriangles; i++)
            {
             
                 centaurFile.read(reinterpret_cast<char*>(&triGlobalNodeIndexes[0]), sizeof(triGlobalNodeIndexes));
                for (int j=0;j<3;j++){globalNodeIndexes[j]=triGlobalNodeIndexes[j]-1;}
                
                addElement(globalNodeIndexes,referenceTriangle);
                
            }
            
            
        }
        centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
        if (checkInt!=sizeOfLine) {std::cerr << "Error in centaur file " << std::endl; return;}
        
    
        
        //Now check the number of quaduratiles in the file
        centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
        int numberOfQuads;			
	    centaurFile.read(reinterpret_cast<char*>(&numberOfQuads),sizeof(numberOfQuads));
        std::cout << "File contains " <<numberOfQuads<<" quaduratile(s)" <<std::endl;
        if (checkInt!=sizeOfLine) {std::cerr << "Error in centaur file " << std::endl; return;}
        
        //Now read the quaduritles in
        centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
        ///\bug Code works but we do not understand why the following integer is throwen away
        centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
        if (numberOfQuads>0)
        {
            unsigned int quadGlobalNodeIndexes[4];
            unsigned int temp;
            std::vector<unsigned int> globalNodeIndexes(4);
            Geometry::ReferenceGeometry<2>* referenceSquare = &Geometry::ReferenceSquare::Instance();
            for (int i=0; i<numberOfQuads; i++)
                {
                    //Reading the vertex indices of each quadrilateral.
                    centaurFile.read(reinterpret_cast<char*>(&quadGlobalNodeIndexes[0]), sizeof(quadGlobalNodeIndexes));
            
                    // renumbering of the vertices to match the ordering assumed by
                    // hpGem:
                    
                    temp=quadGlobalNodeIndexes[2];
                    quadGlobalNodeIndexes[2]=quadGlobalNodeIndexes[3];
                    quadGlobalNodeIndexes[3]=temp;
            
                    // renumber them from 1..N to 0..N-1.
                    for (int j=0; j < 4; j++)
                        globalNodeIndexes[j] = quadGlobalNodeIndexes[j]-1; 
      
                    addElement(globalNodeIndexes,referenceSquare);
                }

	    }
        //Check the line was read correctly : each line in centaur start and end with the line size as a check
        centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
        if (checkInt!=sizeOfLine) {std::cerr << "Error in centaur file " << std::endl; return;}
        
     
        faceFactory();
        
        
    }
    else 
    {
        std::cerr << "Incorrect mesh file. This mesh appears to contain three DIMensional data" << std::endl;
        
    }

    
    
}




/// \bug does not do the bc flags yet or periodic mesh reads
template<unsigned int DIM>
void MeshManipulator<DIM>::faceFactory()
{

    //This will store the number of faces.
    unsigned int numOfFaces;
    
    //Half face holds global node number of both nodes, the element number and the local face number in that element.
    halfFaceDescriptionT halfFace;
    
    //List to hold the half faces. List is used for the quick sorting of the halfFaces
    std::list<halfFaceDescriptionT> halfFaceList;
    
    VectorOfPointIndexesT globalFaceIndexes;
    
    std::vector<Base::Element<DIM>* > tempElementVector(elements_.size());
    
    //first loop over the elements reading in all the data of the half faces
    unsigned int elementID=0;
    for (typename ListOfElementsT::iterator it=elements_.begin(); it !=elements_.end(); ++it)
    {
        const Geometry::PhysicalGeometry<DIM>* const myPhysicalGeometry = (*it)->getPhysicalGeometry();
        const Geometry::ReferenceGeometry<DIM>* const myReferenceGeometry = (*it)->getReferenceGeometry();
        tempElementVector[elementID]=&(**it);
        
        numOfFaces = myReferenceGeometry->getNrOfCodim1Entities();
        
        //Loop over the faces on the element
        for (int i=0; i<numOfFaces; i++)
        {
         
            myPhysicalGeometry->getGlobalFaceNodeIndices(i,globalFaceIndexes);
           // cout<<elementID  <<" " <<globalFaceIndexes[0] <<","<<globalFaceIndexes[1]<< endl;
            
            halfFace.elementNum=elementID;
            halfFace.localFaceIndex=i;
            //Make sure the firstNode number is the smaller of the two global node number, required for the sorting.
            if (globalFaceIndexes[0]<globalFaceIndexes[1])
            {
                halfFace.firstNode=globalFaceIndexes[0];
                halfFace.secondNode=globalFaceIndexes[1];
            }
            else 
            {
                halfFace.firstNode=globalFaceIndexes[1];
                halfFace.secondNode=globalFaceIndexes[0];
            }
            
            //Add the halfFace to the list
            halfFaceList.push_front(halfFace);

            
           
            
        }
            
        elementID++;
  
                 
    }
    
    //Now sort the list on the two value (firstNode, secondNode) so it order and the two pair internal halfFaces are next to each other
    halfFaceList.sort(compareHalfFace);
    
    //Writing out for testing
    for (typename std::list<halfFaceDescriptionT>::const_iterator cit=halfFaceList.begin(); cit !=halfFaceList.end(); ++cit)
    {
        std::cout << (*cit).elementNum << " " << (*cit).localFaceIndex<< " "<< (*cit).firstNode<<" " << (*cit).secondNode<<" "<<  std::endl;
    }
    
    //Now create the faces
    halfFaceDescriptionT current;
    halfFaceDescriptionT next;
    
    //Loop over all the half faces
     for (typename std::list<halfFaceDescriptionT>::const_iterator cit=halfFaceList.begin(); cit !=halfFaceList.end(); ++cit)
     {
         //For the first halfFace just read it in, as we cannot tell if is an internal or external yet.
         if (cit==halfFaceList.begin())
         {
             current=*cit;
         }
         else 
         {
             next=*cit;
             
             //This is an interal face
             if ((current.firstNode==next.firstNode) && (current.secondNode==next.secondNode)) 
             {
                 
                 addFace(tempElementVector[current.elementNum],current.localFaceIndex,tempElementVector[next.elementNum],next.localFaceIndex);
                 //jump a face
                 ++cit;
                 current=*cit;
                 
                 //There was only one face left, so it now the end of the list
                 if (cit==halfFaceList.end())
                 {
                     break;
                 }
                 
                 
             }
             else //it is a boundary face 
             {
                 addFace(tempElementVector[current.elementNum],current.localFaceIndex,NULL,0);
                 current=next;
             }
             
             

         }

         
         
     }
}


  //---------------------------------------------------------------------
template<unsigned int DIM>
int MeshManipulator<DIM>::size() const 
{ 
    return vecOfElementTree_.size();
}

//! Create a new (empty) mesh-tree.
template<unsigned int DIM>
void MeshManipulator<DIM>::createNewMeshTree()
{
  vecOfElementTree_.push_back(new ElementLevelTreeT);
  vecOfFaceTree_.push_back(new FaceLevelTreeT);
  ++numMeshTree_;
  setActiveMeshTree(numMeshTree_ - 1);
}

//! Get the element container of a specific mesh-tree.
template<unsigned int DIM>
typename MeshManipulator<DIM>::ElementLevelTreeT* MeshManipulator<DIM>::ElCont(int meshTreeIdx = -1) const
{
    int onIndex;
    if ((meshTreeIdx >= 0) && (meshTreeIdx < numMeshTree_))
    {
      onIndex = meshTreeIdx;
    }
    else if ((meshTreeIdx == -1) && (activeMeshTree_ >= 0))
    {
      onIndex = activeMeshTree_;
    }
    else
      throw "MeshManipulator::ElCont(): invalid mesh-tree index or no active mesh-tree.";
    
    return vecOfElementTree_[onIndex];
}

//! Get the face container of a specific mesh-tree.
template<unsigned int DIM>
typename MeshManipulator<DIM>::FaceLevelTreeT* MeshManipulator<DIM>::FaCont(int meshTreeIdx = -1) const
{
    int onIndex;
    if ((meshTreeIdx >= 0) && (meshTreeIdx < numMeshTree_))
    {
      onIndex = meshTreeIdx;
    }
    else if ((meshTreeIdx == -1) && (activeMeshTree_ >= 0))
    {
      onIndex = activeMeshTree_;
    }
    else
      throw "MeshManipulator::FaCont(): invalid mesh-tree index or no active mesh-tree.";
    
    return vecOfFaceTree_[onIndex];
}

//! Some mesh generator: centaur / rectangular / triangle / tetrahedra / triangular-prism.
template<unsigned int DIM>
void MeshManipulator<DIM>::someMeshGenerator(int meshTreeIdx = -1)
{
    int onIndex;
    if ((meshTreeIdx >= 0) && (meshTreeIdx < numMeshTree_))
    {
      onIndex = meshTreeIdx;
    }
    else if ((meshTreeIdx == -1) && (activeMeshTree_ >= 0))
    {
      onIndex = activeMeshTree_;
    }
    else
      throw "MeshManipulator::someMeshGenerator(): invalid mesh-tree index or no active mesh-tree.";
      
    const int numberOfElement = 1 + (rand() % 10);
    const int startElementId = (onIndex+1)*1000;
    for (int id=startElementId; id<startElementId+numberOfElement; ++id)
    {
      ElementT el(id);
      vecOfElementTree_[onIndex]->addEntry(el);
    }

    const int numberOfFace = 1 + (rand() % 10);
    const int startFaceId = (onIndex+1)*1000;
    for (int id=startFaceId; id<startFaceId+numberOfFace; ++id)
    {
      FaceT fa(id);
      vecOfFaceTree_[onIndex]->addEntry(fa);
    }
}

//! Set active mesh-tree.
template<unsigned int DIM>
void MeshManipulator<DIM>::setActiveMeshTree(unsigned int meshTreeIdx)
{
    if (meshTreeIdx < numMeshTree_)
        activeMeshTree_ = meshTreeIdx;
    else
        throw "MeshManipulator<DIM>::setActiveMeshTree(): invalid mesh-tree index.\n";
}

//! Get active mesh-tree index.
template<unsigned int DIM>
int MeshManipulator<DIM>::getActiveMeshTree() const
{
    return activeMeshTree_;
}

//! Reset active mesh-tree.
template<unsigned int DIM>
void MeshManipulator<DIM>::resetActiveMeshTree()
{
    activeMeshTree_ = -1;
}

//! Get maximum h-level of a specific mesh-tree.
template<unsigned int DIM>
unsigned int MeshManipulator<DIM>::getMaxLevel(int meshTreeIdx = -1) const
{
    int onIndex;
    if ((meshTreeIdx >= 0) && (meshTreeIdx < numMeshTree_))
    {
      onIndex = meshTreeIdx;
    }
    else if ((meshTreeIdx == -1) && (activeMeshTree_ >= 0))
    {
      onIndex = activeMeshTree_;
    }
    else
      throw "MeshManipulator::getMaxLevel(): invalid mesh-tree index or no active mesh-tree.";
    
    return vecOfElementTree_[onIndex].getMaxLevel();
}

//! Set active level of a specific mesh-tree.
template<unsigned int DIM>
void MeshManipulator<DIM>::setActiveLevel(unsigned int meshTreeIdx, int level)
{
    int onIndex;
    if ((meshTreeIdx >= 0) && (meshTreeIdx < numMeshTree_))
    {
      onIndex = meshTreeIdx;
    }
    else if ((meshTreeIdx == -1) && (activeMeshTree_ >= 0))
    {
      onIndex = activeMeshTree_;
    }
    else
      throw "MeshManipulator::setActiveLevel(): invalid mesh-tree index or no active mesh-tree.";

    if ((level >= 0) && (level <= vecOfElementTree_[onIndex].getMaxLevel()))
    {
      vecOfElementTree_[onIndex].setActiveLevel(level);
      vecOfFaceTree_[onIndex].setActiveLevel(level);
    }
    else
      throw "MeshManipulator::setActiveLevel(): invalid level.";
    
}

//! Get active level of a specific mesh-tree.
template<unsigned int DIM>
int MeshManipulator<DIM>::getActiveLevel(int meshTreeIdx = -1) const
{
    int onIndex;
    if ((meshTreeIdx >= 0) && (meshTreeIdx < numMeshTree_))
    {
      onIndex = meshTreeIdx;
    }
    else if ((meshTreeIdx == -1) && (activeMeshTree_ >= 0))
    {
      onIndex = activeMeshTree_;
    }
    else
      throw "MeshManipulator::getActiveLevel(): invalid mesh-tree index or no active mesh-tree.";

    return vecOfElementTree_[onIndex].getActiveLevel();
}

//! Reset active level of a specific mesh-tree.
template<unsigned int DIM>
void MeshManipulator<DIM>::resetActiveLevel(int meshTreeIdx = -1)
{
    int onIndex;
    if ((meshTreeIdx >= 0) && (meshTreeIdx < numMeshTree_))
    {
      onIndex = meshTreeIdx;
    }
    else if ((meshTreeIdx == -1) && (activeMeshTree_ >= 0))
    {
      onIndex = activeMeshTree_;
    }
    else
      throw "MeshManipulator::resetActiveLevel(): invalid mesh-tree index or no active mesh-tree.";

    vecOfElementTree_[onIndex].resetActiveLevel();
    vecOfFaceTree_[onIndex].resetActiveLevel();
}

//! Duplicate mesh contents including all refined meshes.
template<unsigned int DIM>
void MeshManipulator<DIM>::duplicate(unsigned int fromMeshTreeIdx, unsigned int toMeshTreeIdx, unsigned int upToLevel = 0)
{
    
}

//! Refine a specific mesh-tree.
template<unsigned int DIM>
void MeshManipulator<DIM>::doRefinement(unsigned int meshTreeIdx, int refinementType = -1)
{
    int level = getMaxLevel(meshTreeIdx);
    setActiveLevel(meshTreeIdx,level);

    if (refinementType == -1)  
    {
        // elements should have been flagged before calling this
        std::cout << "MeshManipulator<" << DIM << ">::doRefinement(Mesh(" << meshTreeIdx << "))\n";
        doElementRefinement(meshTreeIdx);
        doFaceRefinement(meshTreeIdx);
        for (ElementIteratorT it=ElCont(meshTreeIdx).beginLevel(level); it != ElCont(meshTreeIdx).end(); ++it)
        {
            // reset element's marking for refinement
//           it->unsetRefineType();
//           it->setBeingRefinedOff();
        }
    }
    else
    {
        // apply the uniform mesh refinement
        for (ElementIteratorT it=ElCont(meshTreeIdx).beginLevel(level); it != ElCont(meshTreeIdx).end(); ++it)
        {
            // mark element for refinement
            it->setRefineType(refinementType);
        }

        std::cout << "MeshManipulator<" << DIM << ">::doRefinement(" << meshTreeIdx << "," << refinementType << ")\n";
        doElementRefinement(meshTreeIdx);
        doFaceRefinement(meshTreeIdx);
        for (ElementIteratorT it=ElCont(meshTreeIdx).beginLevel(level); it != ElCont(meshTreeIdx).end(); ++it)
        {
            // reset element's marking for refinement
//           it->unsetRefineType();
//           it->setBeingRefinedOff();
        }
    }
    
    level = getMaxLevel(meshTreeIdx);
    setActiveLevel(meshTreeIdx,level);
}

//! Do refinement on the elements.
template<unsigned int DIM>
void MeshManipulator<DIM>::doElementRefinement(unsigned int meshTreeIdx)
{
    unsigned int needDummyFaceOnLevel = 0;
    std::vector<ElementT*> vecElementsToRefined;  // list of unrefined elements
    for (ElementIteratorT el=ElCont(meshTreeIdx).begin(); el != ElCont(meshTreeIdx).end(); ++el)
    {
        Geometry::RefinementGeometry<DIM>* RG = el->getRefinementGeometry();
        
        int refineType;
//         refineType = el->getRefineType();   // TODO: add this to Element

/*        

        if (refineType  < 0)   // not refined, just duplicate the element
        {
            vecElementsToRefined.push_back(&(*el));
            needDummyFaceOnLevel = el->getLevel()+1;
            continue;
        }

        //********* Add new nodes to the container
        //----------------------------------------
        RG->getAllNodes(refineType, VectorOfPointPhysicalsT& nodes);
        unsigned int nrNewNodes = nrOfNewNodes(refineType);
        unsigned int nrAllNodes = nodes.size();
        unsigned int nrOldNodes = nrAllNodes - nrNewNodes;
        
        // get physical nodes of this elements
        VectorOfPhysicalPointsT nodesPhys;    // vector of physical points
        VectorOfPointIndexesT   nodesIdx;     // vector of global indices
        for (PointIndexT j=0; j<nrVertices; ++j)
        {
            PointPhysicalT p;       // a physical node
            PointIndexT pIdx;       // a global index
            PG->getVertexPoint (j, p);
            pIdx = PG->getVertexIndex (j);
            nodesPhys.push_back(p);
            nodesIdx.push_back(pIdx);
        }

        // get new physical nodes due to the refinement, and add them
        // up to the nodes collection of this element
        PG->NewPhysicalNodes(refineType, nodesPhys);

        // add the new physical nodes into the nodes container
        // and get their global indices
        for (LocalPointIndexT j=nrVertices; j<nodesPhys.size(); ++j)
        {
            int pIdx = PCptr->getGlobalIndex(nodesPhys[j]);
            if (pIdx >= 0)
            {
                // it's already exist
                nodesIdx.push_back(pIdx);
            }
            else
            {
                // it's not exist yet.  Add it to the node container
                MeshBase<DIM>::AllNodes.addRoot(nodesPhys[j]);
                nodesIdx.push_back(PCptr->size()-1);
            }
        }
        // Now we already have all nodes: being used by this element and
        // to be used by the new sub-elements.

        //********* Add sub-elements to the container
        //----------------------------------------
        unsigned int nrNewElements = PG->nrOfSubElements(refineType);
        std::vector<ElementBase*> vecSubElements;
        for (unsigned int j=0; j<nrNewElements; ++j)
        {
            VectorOfPointIndexesT localNodeIndexes;
            PG->SubElementNodeIndexes(refineType, j, localNodeIndexes);
            ElementDescriptor elDescriptor(localNodeIndexes.size());
            for (VectorOfPointIndexesT k=0; k<localNodeIndexes.size(); ++k)
            {
                elDescriptor.addNode(nodesIdx[localNodeIndexes[k]]);
            }

//             ElementFactory<DIM> elFactory(*this);
//             ElementT *elem = elFactory.makeElement(&elDescriptor);
//             elem->setRefinementType(refineType);
//             vecSubElements.push_back(elem);
        }
        
        // Add the sub-elements as the children of this element
        ElementIteratorT elIt = ElCont(meshTreeIdx).addChildren(el, vecSubElements);
        
        // Clear up the storage
        nodesPhys.clear(); // clear vector of physical points
        nodesIdx.clear();  // clear vector of global index

        
        //********* Add sub-Internal Faces to the container
        //----------------------------------------
        VectorOfPointIndexesT elementIdx1;
        VectorOfPointIndexesT elementIdx2;
        VectorOfPointIndexesT localFaceIdx1;
        VectorOfPointIndexesT localFaceIdx2;
        PG->AdjacentSubElementsPairs(refineType, elementIdx1,localFaceIdx1, elementIdx2,localFaceIdx2);
        unsigned int nrOfNewFaces = elementIdx1.size();

        std::vector<FaceT*> vecSubFaces;
        for (unsigned int j=0; j<nrOfNewFaces; ++j)
        {
            // Create the face, connecting the two elements
            ElementT* el1 = vecSubElements[elementIdx1[j]];
            ElementT* el2 = vecSubElements[elementIdx2[j]];
            FaceT* iFace = new FaceT( el1, localFaceIdx1[j], el2, localFaceIdx2[j]);
            vecSubFaces.push_back(iFace);
        }
        // Add them to the container as the children of a dummy face
        FaceIteratorT fa = FaCont(meshTreeIdx).getDummyFace(el->getLevel());

        FaCont(meshTreeIdx).addChildren(fa, vecSubFaces);
        vecSubFaces.clear();
        
        // Mark that this element was just refined
        el->setJustRefined();
        vecSubElements.clear();
*/
    } // end of loop over elements

    if (needDummyFaceOnLevel)
    {
        FaCont(meshTreeIdx).getDummyFace(needDummyFaceOnLevel);
    }
    
    std::vector<ElementT*> vecEmpty;  // empty list of new elements
    while(!vecElementsToRefined.empty()) 
    {
        ElementT *elem = vecElementsToRefined.back();
        ElementIteratorT elIt = ElCont(meshTreeIdx).addChildren(elem->getIterator(), vecEmpty);
        vecElementsToRefined.pop_back();
    }
}

//! Do refinement on the faces.
template<unsigned int DIM>
void MeshManipulator<DIM>::doFaceRefinement(unsigned int meshTreeIdx)
{
    
}

//! Check whether the two elements may be connected by a face or not.
template<unsigned int DIM>
void MeshManipulator<DIM>::pairingCheck(
                  const ElementIteratorT elL, unsigned int locFaceNrL,
                  const ElementIteratorT elR, unsigned int locFaceNrR,
                  int& pairingValue, bool& sizeOrder)
// pairingValue: 0=not match, 1=partial match, 2=perfect match
// sizeOrder: true = LR,  false = RL
{
    // get node numbers from left (and right) sides
    std::vector<PointIndexT> globNodesL;
    std::vector<PointIndexT> globNodesR;
    
    const Geometry::PhysicalGeometry<DIM>* const leftPG = elL->getPhysicalGeometry();
    Geometry::PhysicalGeometry<DIM>* rightPG;
    
    leftPG->getFaceNodeIndices(locFaceNrL,globNodesL);
    if (&*elR != 0)
    {
        rightPG = elR->getPhysicalGeometry();
        rightPG->getFaceNodeIndices(locFaceNrR,globNodesR);
    }

    // store them as sets
    std::set<PointIndexT> setL(globNodesL.data(), globNodesL.data() + globNodesL.size());
    std::set<PointIndexT> setR(globNodesR.data(), globNodesR.data() + globNodesR.size());

    if (setL == setR)
    {
        // nodes of the faces are match
        pairingValue = 2;
        sizeOrder = true;
        return;
    }
    
    std::set<PointIndexT> commNodes;
    std::set_intersection( setL.begin(), setL.end(), setR.begin(), setR.end(), std::inserter(commNodes, commNodes.begin()) );

    unsigned int nrNodesL = globNodesL.size();
    unsigned int nrNodesR = globNodesR.size();
    if (nrNodesL != nrNodesR)
    {
        pairingValue = 0;
        return;
    }
    
    std::set<PointIndexT> diffNodesL;
    std::set_difference( setL.begin(), setL.end(), commNodes.begin(), commNodes.end(), std::inserter(diffNodesL, diffNodesL.end()) );

    std::set<PointIndexT> diffNodesR;
    std::set_difference( setR.begin(), setR.end(), commNodes.begin(), commNodes.end(), std::inserter(diffNodesR, diffNodesR.end()) );

    typedef std::set<PointIndexT>::iterator SetIterType;

    // do collinear test for all possible pairs
    pairingValue = 1;
    
    // order of face sizes: true = LR;  false = RL
    sizeOrder = true;
    
    // loop until empty or until the faces proved not collinear
    while (!diffNodesL.empty())
    {
        SetIterType itDiffL=diffNodesL.begin();
        PointPhysicalT ppL;
        leftPG->getGlobalNodeCoordinates(*itDiffL, ppL);

        bool foundCollinear(false);
        SetIterType itDiffL2;
        SetIterType itDiffR;
        SetIterType itDiffR2;
        for (itDiffR=diffNodesR.begin(); itDiffR!=diffNodesR.end(); ++itDiffR)
        {
          PointPhysicalT ppR;
          rightPG->getGlobalNodeCoordinates(*itDiffR, ppR);

          if (commNodes.size() > 0)
          {
            for (SetIterType itCommon=commNodes.begin(); itCommon!=commNodes.end(); ++itCommon)
            {
              PointPhysicalT pp0;
              leftPG->getGlobalNodeCoordinates(*itCommon, pp0);

              //-------------
              // perform 3-nodes collinear test
              PointPhysicalT dL(ppL-pp0);
              PointPhysicalT dR(ppR-pp0);

              double ratio = 0.;
              unsigned int d;
              for (d=0; d<DIM; ++d)
              {
                if (std::abs(dR[d])>Geometry::SmallerDoubleThanMinimalSizeOfTheMesh)
                {
                  ratio = dL[d]/dR[d];
                  break;
                }
              }

              if (ratio < 0.)
              {
                pairingValue = 0;
                return;
              }
              
              foundCollinear = true;
              for (;d<DIM; ++d)
              {
                if (std::abs(dR[d])>Geometry::SmallerDoubleThanMinimalSizeOfTheMesh)
                {
                  if (dL[d]/dR[d] < 0.)
                  {
                    foundCollinear = false;
                    break;
                  }

                  if (std::abs(dL[d]/dR[d] - ratio) > Geometry::SmallerDoubleThanMinimalSizeOfTheMesh)
                  {
                    foundCollinear = false;
                    break;
                  }
                }
              }

              // order of face sizes: true = LR;  false = RL
              sizeOrder = (sizeOrder && (ratio <= 1.0));
              
              // found a collinear nodes pair
              if (foundCollinear) break;
            }  // for itCommon
          } // if (commNodes.size() > 0)
          else
          // no common nodes
          {
            for (itDiffL2=diffNodesL.begin(); itDiffL2!=diffNodesL.end(); ++itDiffL2)
            {
              if ((itDiffL2 == itDiffL) || (itDiffL2 == itDiffR)) 
                  continue;
              
              PointPhysicalT ppL2;
              leftPG->getGlobalNodeCoordinates(*itDiffL2, ppL2);

              for (itDiffR2=diffNodesR.begin(); itDiffR2!=diffNodesR.end(); ++itDiffR2)
              {
                if ((*itDiffR2 == *itDiffR) || (*itDiffR2 == *itDiffL) || (*itDiffR2 == *itDiffL2)) 
                    continue;
                
                PointPhysicalT ppR2;
                rightPG->getGlobalNodeCoordinates(*itDiffR2, ppR2);
              
                //-------------
                // perform two 3-nodes pairs collinear test 
                PointPhysicalT dL1(ppL-ppR);
                PointPhysicalT dR1(ppL2-ppR);
                
                PointPhysicalT dL2(ppL2-ppR);
                PointPhysicalT dR2(ppR2-ppR);
                
                PointPhysicalT dL3(ppL-ppR2);
                PointPhysicalT dR3(ppR-ppR2);
                
                double ratio1 = 0.;
                double ratio2 = 0.;
                double ratio3 = 0.;
                unsigned int d1;
                unsigned int d2;
                unsigned int d3;
                for (d1=0; d1<DIM; ++d1)
                {
                  if (std::abs(dR1[d1])>Geometry::SmallerDoubleThanMinimalSizeOfTheMesh)
                  {
                    ratio1 = dL1[d1]/dR1[d1];
                    break;
                  }
                }

                for (d2=0; d2<DIM; ++d2)
                {
                  if (std::abs(dR2[d2])>Geometry::SmallerDoubleThanMinimalSizeOfTheMesh)
                  {
                    ratio2 = dL2[d2]/dR2[d2];
                    break;
                  }
                }
                
                for (d3=0; d3<DIM; ++d3)
                {
                  if (std::abs(dR3[d3])>Geometry::SmallerDoubleThanMinimalSizeOfTheMesh)
                  {
                    ratio3 = dL3[d3]/dR3[d3];
                    break;
                  }
                }
                
                if ((ratio1 < 0) || (ratio2 < 0) || (ratio3 < 0))
                {
                  pairingValue = 0;
                  return;
                }

                foundCollinear = true;
                for (;d1<DIM; ++d1)
                {
                  if (std::abs(dR1[d1])>Geometry::SmallerDoubleThanMinimalSizeOfTheMesh)
                  {
                    if (dL1[d1]/dR1[d1] < 0.)
                    {
                      pairingValue = 0;
                      return;
                    }

                    if (std::abs(dL1[d1]/dR1[d1] - ratio1) > Geometry::SmallerDoubleThanMinimalSizeOfTheMesh)
                    {
                      foundCollinear = false;
                      break;
                    }
                  }
                }
                
                for (;(d2<DIM) && foundCollinear; ++d2)
                {
                  if (std::abs(dR2[d2])>Geometry::SmallerDoubleThanMinimalSizeOfTheMesh)
                  {
                    if (dL2[d2]/dR2[d2] < 0.)
                    {
                      pairingValue = 0;
                      return;
                    }

                    if (std::abs(dL2[d2]/dR2[d2] - ratio2) > Geometry::SmallerDoubleThanMinimalSizeOfTheMesh)
                    {
                      foundCollinear = false;
                      break;
                    }
                  }
                }
                
                for (;(d3<DIM) && foundCollinear; ++d3)
                {
                  if (std::abs(dR3[d3])>Geometry::SmallerDoubleThanMinimalSizeOfTheMesh)
                  {
                    if (dL3[d3]/dR3[d3] < 0.)
                    {
                      pairingValue = 0;
                      return;
                    }

                    if (std::abs(dL3[d3]/dR3[d3] - ratio3) > Geometry::SmallerDoubleThanMinimalSizeOfTheMesh)
                    {
                      foundCollinear = false;
                      break;
                    }
                  }
                }
                
                // found a collinear nodes pair
                if (foundCollinear) break;     // itDiffR2 will be deleted
              }
              // found a collinear nodes pair
              if (foundCollinear) break;     // itDiffL2 will be deleted
            }  // for itDiffL2
          }  //  else, no common nodes
         
          // found a collinear nodes pair
          if (foundCollinear) break;  // itDiffR will be deleted
        }  // for itDiffR

        if (foundCollinear) 
        // found a collinear nodes pair
        {
          if (commNodes.size() > 0)
          {
            diffNodesR.erase(itDiffR);
          }
          else
          {
            diffNodesL.erase(itDiffL2);
            diffNodesR.erase(itDiffR);
            diffNodesR.erase(itDiffR2);
          }
        }
        else
        // the faces are not collinear
        {
          pairingValue = 0;
          break;
        }
        
        diffNodesL.erase(itDiffL);
    }  // end while-loop
    
    if (!diffNodesR.empty())
      pairingValue = 0;

}  // end pairingCheck

//! Check whether the two elements may be connected by a face or not in periodic face case.
template<unsigned int DIM>
void MeshManipulator<DIM>::periodicPairingCheck(const FaceIteratorT fa, 
                  const ElementIteratorT elL, unsigned int localFaceNrL,
                  const ElementIteratorT elR, unsigned int localFaceNrR,
                  int& pairingValue, bool& sizeOrder)
{
    
}
  //---------------------------------------------------------------------
