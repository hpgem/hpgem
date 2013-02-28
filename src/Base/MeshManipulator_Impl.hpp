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

    

    
    //This stores the number of nodes in each codimension i.e. if you have 2 by 2 element it is 3 nodes 
    std::vector<int> numOfNodesInEachSubspace(DIM), numOfElementsInEachSubspace(DIM);
    
    numOfNodesInEachSubspace[0]=1;
    numOfElementsInEachSubspace[0]=1;
    
    //This will be the total number of nodes required in the problem
    unsigned int totalNumOfNodes,totalNumOfElements, verticesPerElement;
    
    totalNumOfNodes=(linearNoElements[0]+1);
    totalNumOfElements=(linearNoElements[0]);
    verticesPerElement=2;
    int powerOf2;
    
    for (int idim=1; idim<DIM;++idim)
    {
        totalNumOfNodes*=(linearNoElements[idim]+1);
        totalNumOfElements*=(linearNoElements[idim]);
        verticesPerElement*=2;
        
        numOfElementsInEachSubspace[idim]=numOfElementsInEachSubspace[idim-1]*(linearNoElements[idim-1]);
        numOfNodesInEachSubspace[idim]=numOfNodesInEachSubspace[idim-1]*(linearNoElements[idim-1]+1);
        
    }
    
    //temp point for storing the node locations
    PointPhysicalT x;
    
    
    ///Stage 2 : Create the nodes
    //Now loop over all the nodes and calculate the coordinates for rach dimension (this makes the algorithm independent of dimension
    for (int nodeIndex=0; nodeIndex<totalNumOfNodes; ++nodeIndex)
    {
        int nodeIndexRemain=nodeIndex;
        
        
        
        for (int idim=DIM-1;idim>-1;--idim)
        {
            x[idim] = BottomLeft[idim] + (nodeIndexRemain/numOfNodesInEachSubspace[idim]*delta_x[idim]);
            nodeIndexRemain %=numOfNodesInEachSubspace[idim];
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
    
    //This will store a two array, which goes elementNumber and the n-dimensional coordinate of that element. This is only used for the face creation (step 4).
    //std::vector<std::vector<unsigned int> > allElementNdId;
    //allElementNdId.resize(totalNumOfElements);
    
  
    //elementNdId is DIM coordinate of the bottom left node i.e. in two (0,0), (1,0) ,(2,0) ... etc are the first three (if at least three elements in x)
    for (int elementIndex=0; elementIndex<totalNumOfElements; ++elementIndex)
    {
        int numElementsRemaining=elementIndex;
        
        for (int idim=DIM-1; idim>-1;--idim)
        {
            elementNdId[idim]=numElementsRemaining/numOfElementsInEachSubspace[idim];
            numElementsRemaining %= numOfElementsInEachSubspace[idim];
        }
    //    allElementNdId[elementIndex].resize(DIM);
    //    allElementNdId[elementIndex]=elementNdId;
        
        // vertexNdId are the DIM coordinate of each vertex in the element with vertexNdId[0] being the bottom left
        for (int i=0; i<verticesPerElement; ++i)
        {
            powerOf2 = 1;
            for (int idim=0; idim<DIM; ++idim)
            {
                vertexNdId[idim] = elementNdId[idim] + ((i & powerOf2) !=0);
                powerOf2 *=2;
            }
            globalVertexID[i]=vertexNdId[0];
                
            //Now map to the one dimensional global ID
            for (int idim=1;idim<DIM;++idim){globalVertexID[i] += vertexNdId[idim]*numOfNodesInEachSubspace[idim];}
        }
            
        tempElementVector[elementIndex]=addElement(globalVertexID,referenceGeometry);
    }
    
    
    //Stage 4 : Create the faces
    /// \bug Only implemented for 1D and 2D at the momennt
    // Note the linear number of faces in each direction is linearNoElements-1;
    
    //loop over all the elements
    //for (int elementIndex=0;elementIndex<totalNumOfElements;elementIndex++)
//    {
//        //loop over the dimensions of the element. Element create internal faces in possitive direction and boundary faces
//        for (int dimIndex=0; dimIndex<DIM;dimIndex++)
//        {
//            //Left x boundary
//            if (allElementNdId[elementIndex][dimIndex]==0)
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
       
     //   default:
     //       throw("Face generator not implemented in this dimension");
            
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
            std::cerr<<"Centaur mesh reader has not been implemented in this dimension" << std::endl;
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
    
    // Centaur File Type <0 is two dimensional and >0 is three dimensional
    int centaurFileType;			
	centaurFile.read(reinterpret_cast<char*>(&centaurFileType),sizeof(centaurFileType));
    
    
    
    if (centaurFileType<0)
    {
        std::cout << "Reading a two dimensional centaur mesh" <<std::endl;
        
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
        std::cerr << "Incorrect mesh file. This mesh appears to contain three dimensional data" << std::endl;
        
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
typename MeshManipulator<DIM>::ElementLevelTreeT* MeshManipulator<DIM>::ElCont(int meshTreeIdx) const
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
typename MeshManipulator<DIM>::FaceLevelTreeT* MeshManipulator<DIM>::FaCont(int meshTreeIdx) const
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
void MeshManipulator<DIM>::someMeshGenerator(int meshTreeIdx)
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
unsigned int MeshManipulator<DIM>::getMaxLevel(int meshTreeIdx) const
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
int MeshManipulator<DIM>::getActiveLevel(int meshTreeIdx) const
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
void MeshManipulator<DIM>::resetActiveLevel(int meshTreeIdx)
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
void MeshManipulator<DIM>::duplicate(unsigned int fromMeshIdx, unsigned int toMeshIdx, unsigned int upToLevel)
{
}

//! Refine a specific mesh-tree.
template<unsigned int DIM>
void MeshManipulator<DIM>::doRefinement(unsigned int meshTreeIdx, int refinementType)
{
}

//! Do refinement on the elements.
template<unsigned int DIM>
void MeshManipulator<DIM>::doElementRefinement(unsigned int meshTreeIdx)
{
}

//! Do refinement on the faces.
template<unsigned int DIM>
void MeshManipulator<DIM>::doFaceRefinement(unsigned int meshTreeIdx)
{
}

//! Check whether the two elements may be connected by a face or not.
template<unsigned int DIM>
void MeshManipulator<DIM>::pairingCheck(const ElementIteratorT elL, unsigned int locFaceNrL,
                  const ElementIteratorT elR, unsigned int locFaceNrR,
                  int& pairingValue, bool& sizeOrder)
{
}

//! Check whether the two elements may be connected by a face or not in periodic face case.
template<unsigned int DIM>
void MeshManipulator<DIM>::periodicPairingCheck(const FaceIteratorT fa, 
                  const ElementIteratorT elL, unsigned int localFaceNrL,
                  const ElementIteratorT elR, unsigned int localFaceNrR,
                  int& pairingValue, bool& sizeOrder)
{
}
  //---------------------------------------------------------------------
