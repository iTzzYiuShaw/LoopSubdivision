///////////////////////////////////////////////////
//
//  Hamish Carr
//  September, 2020
//
//  ------------------------
//  DirectedEdgeSurface.cpp
//  ------------------------
//  
//  Base code for rendering assignments.
//
//  Minimalist (non-optimised) code for reading and 
//  rendering an object file
//  
//  We will make some hard assumptions about input file
//  quality. We will not check for manifoldness or 
//  normal direction, &c.  And if it doesn't work on 
//  all object files, that's fine.
//
//  While I could set it up to use QImage for textures,
//  I want this code to be reusable without Qt, so I 
//  shall make a hard assumption that textures are in 
//  ASCII PPM and use my own code to read them
//  
///////////////////////////////////////////////////

// include the header file
#include "DirectedEdgeSurface.h"

// include the C++ standard libraries we want
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cmath>
#include <fstream>

// include the Cartesian 3- vector class
#include "Cartesian3.h"
#include "SphereVertices.h"

#define MAXIMUM_LINE_LENGTH 1024

// constructor will initialise to safe values
DirectedEdgeSurface::DirectedEdgeSurface()
    : centreOfGravity(0.0,0.0,0.0)
    { // DirectedEdgeSurface()
    // force arrays to size 0
    vertices.resize(0);
    normals.resize(0);
	firstDirectedEdge.resize(0);
	faceVertices.resize(0);
	otherHalf.resize(0);
    } // DirectedEdgeSurface()

// read routine returns true on success, failure otherwise
bool DirectedEdgeSurface::ReadObjectStream(std::istream &geometryStream)
    { // ReadObjectStream()
    
    // create a read buffer
    char readBuffer[MAXIMUM_LINE_LENGTH];
    
    // the rest of this is a loop reading lines & adding them in appropriate places
    while (true)
        { // not eof
		// token for identifying meaning of line
		std::string token;

        // character to read
        geometryStream >> token;
        
        // check for eof() in case we've run out
        if (geometryStream.eof())
            break;

        // otherwise, switch on the token we read
		if (token == "#")
			{ // comment 
			// read and discard the line
			geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
            } // comment
		else if (token == "Vertex")
			{ // vertex
			// variables for the read
			unsigned int vertexID;
			geometryStream >> vertexID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
			if (vertexID != vertices.size())
				{ // bad vertex ID
				// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
				} // bad vertex ID				
			
			// read in the new vertex position
			Cartesian3 newVertex;
			geometryStream >> newVertex;
			
			// and add it to the vertices
			vertices.push_back(newVertex);
			} // vertex
		else if (token == "Normal")
			{ // normal
			// variables for the read
			unsigned int normalID;
			geometryStream >> normalID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
			if (normalID != normals.size())
				{ // bad ID
				// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
				} // bad ID				
			
			// read in the new normal
			Cartesian3 newNormal;
			geometryStream >> newNormal;
			
			// and add it to the vertices
			normals.push_back(newNormal);
			} // normal
		else if (token == "FirstDirectedEdge")
			{ // first directed edge
			// variables for the read
			unsigned int FDEID;
			geometryStream >> FDEID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
			if (FDEID != firstDirectedEdge.size())
				{ // bad ID
				// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
				} // bad ID				
			
			// read in the new FDE
			unsigned int newFDE;
			geometryStream >> newFDE;
			
			// and add it to the vertices
			firstDirectedEdge.push_back(newFDE);
			//std::cout<< "FH: " << newFDE << std::endl;
			} // first directed edge
		else if (token == "Face")
			{ // face
			// variables for the read
			unsigned int faceID;
			geometryStream >> faceID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
			if (faceID != faceVertices.size()/3)
				{ // bad face ID
				// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
				} // bad face ID				
			
			// read in the new face vertex (3 times)
			unsigned int newFaceVertex;
			
			geometryStream >> newFaceVertex;
			//std::cout <<"face: " << newFaceVertex;
			faceVertices.push_back(newFaceVertex);
			geometryStream >> newFaceVertex;
			//std::cout << newFaceVertex;
			faceVertices.push_back(newFaceVertex);
			geometryStream >> newFaceVertex;
			//std::cout << newFaceVertex<<std::endl;
			faceVertices.push_back(newFaceVertex);

			} // face
		else if (token == "OtherHalf")
			{ // other half
			// variables for the read
			unsigned int otherHalfID;
			geometryStream >> otherHalfID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
			if (otherHalfID != otherHalf.size())
				{ // bad ID
				// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
				} // bad ID				
			
			// read in the new face vertex (3 times)
			unsigned int newOtherHalf;
			geometryStream >> newOtherHalf;
			otherHalf.push_back(newOtherHalf);
			//std::cout <<"Ot: " << newOtherHalf << std::endl;
			} // other half
        } // not eof

    // compute centre of gravity
    // note that very large files may have numerical problems with this
    centreOfGravity = Cartesian3(0.0, 0.0, 0.0);

    // if there are any vertices at all
    if (vertices.size() != 0)
        { // non-empty vertex set
        // sum up all of the vertex positions
        for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
            centreOfGravity = centreOfGravity + vertices[vertex];
        
        // and divide through by the number to get the average position
        // also known as the barycentre
        centreOfGravity = centreOfGravity / vertices.size();

        // start with 0 radius
        objectSize = 0.0;

        // now compute the largest distance from the origin to a vertex
        for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
            { // per vertex
            // compute the distance from the barycentre
            float distance = (vertices[vertex] - centreOfGravity).length();         
            
            // now test for maximality
            if (distance > objectSize)
                objectSize = distance;
            } // per vertex
        } // non-empty vertex set

    // return a success code
    return true;
    } // ReadObjectStream()

// write routine
void DirectedEdgeSurface::WriteObjectStream(std::ostream &geometryStream)
    { // WriteObjectStream()
	geometryStream << "#" << std::endl; 
	geometryStream << "# Created for Leeds COMP 5821M Autumn 2020" << std::endl; 
	geometryStream << "#" << std::endl; 
	geometryStream << "#" << std::endl; 
	geometryStream << "# Surface vertices=" << vertices.size() << " faces=" << faceVertices.size()/3 << std::endl; 
	geometryStream << "#" << std::endl; 

	// output the vertices
    for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
        geometryStream << "Vertex " << vertex << " " << std::fixed << vertices[vertex] << std::endl;

    // and the normal vectors
    for (unsigned int normal = 0; normal < normals.size(); normal++)
        geometryStream << "Normal " << normal << " " << std::fixed << normals[normal] << std::endl;

	// and the first directed edges
    for (unsigned int vertex = 0; vertex < firstDirectedEdge.size(); vertex++)
        geometryStream << "FirstDirectedEdge " << vertex<< " " << std::fixed << firstDirectedEdge[vertex] << std::endl;

    // and the faces - increment is taken care of internally
    for (unsigned int face = 0; face < faceVertices.size(); )
        { // per face
        geometryStream << "Face  " << face << "  ";
        
        // read in three vertices
        geometryStream << faceVertices[face++] << "  ";
        geometryStream << faceVertices[face++] << "  ";
        geometryStream << faceVertices[face++];
            
        geometryStream << std::endl;
        } // per face

	// and the other halves
	for (unsigned int dirEdge = 0; dirEdge < otherHalf.size(); dirEdge++)
		geometryStream << "OtherHalf " << dirEdge << "   " << otherHalf[dirEdge] << std::endl;
    } // WriteObjectStream()

// routine to render
void DirectedEdgeSurface::Render(RenderParameters *renderParameters)
    { // Render()
    // Ideally, we would apply a global transformation to the object, but sadly that breaks down
    // when we want to scale things, as unless we normalise the normal vectors, we end up affecting
    // the illumination.  Known solutions include:
    // 1.   Normalising the normal vectors
    // 2.   Explicitly dividing the normal vectors by the scale to balance
    // 3.   Scaling only the vertex position (slower, but safer)
    // 4.   Not allowing spatial zoom (note: sniper scopes are a modified projection matrix)
    //
    // Inside a game engine, zoom usually doesn't apply. Normalisation of normal vectors is expensive,
    // so we will choose option 2.  

    // Scale defaults to the zoom setting
    float scale = renderParameters->zoomScale;
    scale /= objectSize;
        
    //  now scale everything
    glScalef(scale, scale, scale);

    // apply the translation to the centre of the object if requested
    glTranslatef(-centreOfGravity.x, -centreOfGravity.y, -centreOfGravity.z);

    // start rendering
    glBegin(GL_TRIANGLES);

	// set colour for pick render - ignored for regular render
	glColor3f(1.0, 1.0, 1.0);

    // loop through the faces
	for (unsigned int face = 0; face < faceVertices.size(); face +=3)
		{ // per face
		// if we want flat normals, compute them here
		if (renderParameters->useFlatNormals)
			{ // flat normals
			// find two vectors along edges of the triangle
			Cartesian3 pq = vertices[faceVertices[face+1]] - vertices[faceVertices[face]];
			Cartesian3 pr = vertices[faceVertices[face+2]] - vertices[faceVertices[face]];

			// take their cross product and normalise
			Cartesian3 faceNormal = pq.cross(pr).unit();

			// and use it to set the glNormal
			glNormal3f(faceNormal.x * scale, faceNormal.y * scale, faceNormal.z * scale);
			} // flat normals

		// we have made a HARD assumption that we have enough normals
		for (unsigned int vertex = face; vertex < face+3; vertex++)
			{ // per vertex
		
			// if we are using smooth normals
			if (!renderParameters->useFlatNormals)
				// set the normal vector
				glNormal3f
					(
					normals[faceVertices[vertex]].x * scale,
					normals[faceVertices[vertex]].y * scale,
					normals[faceVertices[vertex]].z * scale
					);
			
			// and set the vertex position
			glVertex3f
				(
				vertices[faceVertices[vertex]].x,
				vertices[faceVertices[vertex]].y,
				vertices[faceVertices[vertex]].z
				);

			} // per vertex

		} // per face

    // close off the triangles
    glEnd();
    
    // now we add a second loop to render the vertices if desired
    if (!renderParameters->showVertices)
    	return;

	glDisable(GL_LIGHTING);

	// loop through the vertices
	for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
		{ // per vertex
		// use modelview matrix (not most efficient solution, but quickest to code)
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glTranslatef(vertices[vertex].x, vertices[vertex].y, vertices[vertex].z);
		glScalef(0.1 * renderParameters->vertexSize, 0.1 * renderParameters->vertexSize, 0.1 * renderParameters->vertexSize);
		renderTriangulatedSphere();
		glPopMatrix();
		} // per vertex 
    
    } // Render()

// subdivide the mesh
void DirectedEdgeSurface::Subdivision()
{
	int initialSize = faceVertices.size();
	int iniMaxVerIndex = vertices.size()-1;
	int numOfNewIndex = 0;
	int numOfProcFace = 0;
	int curMaxVerIndex = faceVertices.size()-1;
	int sizeOfOriVer = iniMaxVerIndex+1;

	//Generating new vertices for each edge
	//Get indices for new vertices
	for (int i = 0; i < initialSize; i+=3)
	{//Face

		//Iterate each edge, getting 3 vertices
		//std::cout <<"start"<<std::endl;
		std::vector<unsigned int> newIndices;
		for (int k = 0; k < 3; k++)
		{//Edge
			
			//Check if otherHalf has already been signed a new vertex
			int otherH = otherHalf[i+k];

			//The otherhalf has already been signed a new vertex
			//std::cout << (int) ((float)otherH / 3.0f) <<std::endl;
			if ((int) ((float)otherH / 3.0f) <= numOfProcFace)
			{
				//std::cout << otherH << std::endl;
				newIndices.push_back(faceVertices[otherH + initialSize]);
			}else
			{
				numOfNewIndex++;
				newIndices.push_back(iniMaxVerIndex + numOfNewIndex);
			}
		
		}//Edge
		//std::cout <<"done"<<std::endl;

		int oriIndex_1 = faceVertices[i];
		int oriIndex_2 = faceVertices[i+1];
		int oriIndex_3 = faceVertices[i+2];
		//Pushing faces at corner into the array

		//First corner
		faceVertices.push_back(newIndices[0]);
		faceVertices.push_back(oriIndex_1);
		faceVertices.push_back(newIndices[1]);

		//Second corner
		faceVertices.push_back(newIndices[1]);
		faceVertices.push_back(oriIndex_2);
		faceVertices.push_back(newIndices[2]);

		//Third corner
		faceVertices.push_back(newIndices[2]);
		faceVertices.push_back(oriIndex_3);
		faceVertices.push_back(newIndices[0]);

		//Insert central face
		//The indices of central face is incerted at the position following the original indices tuple
		//This design keeps the original index tuples, which are used to calculate new coordinates for original points
		//and to calculate the coordinates of new points.
		for (int k = 1; k <= 3; k++)
		{
			faceVertices.insert(faceVertices.begin()+(curMaxVerIndex+k), newIndices[k-1]);
		}

		curMaxVerIndex+=3;
		numOfProcFace++;
	}//Facex

	//Updating the coordinate of every new point
	for (int i = initialSize; i < 2 * initialSize; i++)
	{
		//check if current new point is already calculated
		if (faceVertices[i] > vertices.size()-1)
		{			
			//This new point has not been updated
			int indexOfEdge = i - initialSize;
			int indexOfOtherH = otherHalf[indexOfEdge];

			int startOfEdge = indexOfEdge - (indexOfEdge % 3);
			int startOfOtherH = indexOfOtherH - (indexOfOtherH % 3);

			//Get the vertices that form the edge that this new point is loacated at
			Cartesian3 ver1 = vertices[faceVertices[indexOfEdge]];
			Cartesian3 ver2 = vertices[faceVertices[indexOfOtherH]];

			Cartesian3 norOfVer1 = normals[faceVertices[indexOfEdge]];
			Cartesian3 norOfVer2 = normals[faceVertices[indexOfOtherH]];

			//Get diagonal vertices. The ones that are not equal to ver1 and ver2
			Cartesian3 ver3;
			Cartesian3 ver4;
			
			Cartesian3 norOfVer3;
			Cartesian3 norOfVer4;

			//Looking for diagonal vertices
			int inOfVer4;
			for (int k = 0; k < 3; k++)
			{
				if (faceVertices[startOfOtherH + k] != faceVertices[indexOfEdge]&&
					faceVertices[startOfOtherH + k] != faceVertices[indexOfOtherH]
					)

				{
					ver4 = vertices[faceVertices[startOfOtherH + k]];
					norOfVer3 = normals[faceVertices[startOfOtherH + k]];
					inOfVer4 = startOfOtherH + k;
				}
			}

			for (int k = 0; k < 3; k++)
			{
				if (faceVertices[startOfEdge + k] != faceVertices[indexOfEdge]&&
					faceVertices[startOfEdge + k] != faceVertices[indexOfOtherH]&&
					faceVertices[startOfEdge + k] != faceVertices[inOfVer4]
					)
				{
					ver3 = vertices[faceVertices[startOfEdge + k]];
					norOfVer4 = normals[faceVertices[startOfEdge + k]];
				}
			}
				
			Cartesian3 verOfnewPoint = (3.0f/8.0f) * (ver1 + ver2) + (1.0f/8.0f) * (ver3 + ver4);
			Cartesian3 norOfNewPoint = (3.0f/8.0f) * (norOfVer1 + norOfVer2) + (1.0f/8.0f) * (norOfVer3 + norOfVer4);

			//Every new vertex is pushed into the vertices array, the original vertices are kept
			//This design preserve the coordinates of original vertices, which are used to updating coordinates.
			vertices.push_back(verOfnewPoint);
			normals.push_back(norOfNewPoint);
		}
	}

	//Update the coordinate of every original point
	for (int i = 0; i < iniMaxVerIndex+1; i++)
	{	
		Cartesian3 oriVer = vertices[i];
		Cartesian3 oriNor = normals[i];
		std::vector<Cartesian3> adjacentVer;
		std::vector<Cartesian3> norOfAdjacVer;

		//Looking for triangles that include this vertex
		for (int k = 0; k < initialSize; k++)
		{
			//If there is an edge pointing at this vertex
			if (faceVertices[k] == i)
			{
				//Locating the first vertex of this triangle
				int startOfTri = k - (k%3);

				//Loop through 3 edges, push the adjacent vers
				for(int q = 0; q < 3; q++)
				{
					if(faceVertices[startOfTri+q] != i){
						
						bool isRepeated = false;
						for(int adj = 0 ; adj < adjacentVer.size();adj++)
						{
							if(adjacentVer[adj] == vertices[faceVertices[startOfTri+q]])
							{
								isRepeated = true;
								break;
							}
						}
						if(!isRepeated)
						{
							adjacentVer.push_back(vertices[faceVertices[startOfTri+q]]);
							norOfAdjacVer.push_back(normals[faceVertices[startOfTri+q]]);
							//std::cout << "adjacent: " <<faceVertices[startOfTri+q]<<std::endl;
						}
					}
				}
			}
		}

		float n = (float)adjacentVer.size();
		float PI = 3.1415926;
		float a = (1.0f/n) * (5.0f/8.0f - pow((3.0f/8.0f + (1.0f/4.0f)*cos(2.0f*PI/n)),2.0));

		if(n == 3)
			a = 3.0f/16.0f;

		Cartesian3 accum;
		Cartesian3 accumNor;
		for (int k = 0; k < adjacentVer.size(); k++)
		{
			accum = accum + a * adjacentVer[k];
			accumNor = accumNor + a * norOfAdjacVer[k];
		}
		
		Cartesian3 newOri = (1.0f - n*a)*oriVer + accum;
		Cartesian3 newNor = (1.0f - n*a)*oriNor + accumNor;
		
		//All updated original vertex is inserted at the slot following the original vertices.
		//This design preserve the coordinates of original vertices, which are used to calculate
		vertices.insert(vertices.begin()+sizeOfOriVer+i, newOri);
		normals.insert(normals.begin()+sizeOfOriVer+i, newNor);
	}

	//Pop out original vertices
	for (int i = 0; i < sizeOfOriVer; i++)
	{
		vertices.erase(vertices.begin());
	}
	
	//Pop out original normals
	for (int i = 0; i < sizeOfOriVer; i++)
	{
		normals.erase(normals.begin());
	}

	//Pop out indices of original vertices
	for (int i = 0; i < initialSize; i++)
	{
		faceVertices.erase(faceVertices.begin());
	}
	
	otherHalf.clear();
	otherHalf.resize(faceVertices.size(),vertices.size()+2);
	int identifier = vertices.size()+2;
	
	//Updating otherHalf
	for(int i = 0; i < faceVertices.size(); i++)
	{
		if(otherHalf[i] != identifier)
			continue;

		int currentId = faceVertices[i];
		int formerId;
		if(i % 3 == 0)
		{
			formerId = faceVertices[i + 2];
		}else{
			formerId = faceVertices[i - 1];
		}

		for(int k = 0; k < faceVertices.size(); k++)
		{
			
			int tempCur = faceVertices[k];
			int tempFormer;
			if(k % 3 == 0)
			{
				tempFormer = faceVertices[k + 2];
			}else{
				tempFormer = faceVertices[k - 1];
			}
		
			if((formerId == tempCur) && (currentId == tempFormer))
			{
				if (otherHalf[k] == identifier)
				{
					otherHalf[i] = k;
					otherHalf[k] = i;		
				}
			}
		}
	}
}

// output a new directed edge data structure in a new text file
void DirectedEdgeSurface::WriteDirFile(char* fileName)
{
	std::ofstream outFile(fileName,std::ios::out);

	outFile << "#" << std::endl; 
	outFile << "# Created for Leeds COMP 5821M Autumn 2020" << std::endl; 
	outFile << "#" << std::endl; 
	outFile << "#" << std::endl; 
	outFile << "# Surface vertices=" << vertices.size() << " faces=" << faceVertices.size()/3 << std::endl; 
	outFile << "#" << std::endl; 

	// write the vertices
    for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
        outFile << "Vertex " << vertex << "  " << std::fixed << vertices[vertex] << std::endl;

    // and the normal vectors
    for (unsigned int normal = 0; normal < normals.size(); normal++)
        outFile << "Normal " << normal << "  " << std::fixed << normals[normal] << std::endl;

	// and the first directed edges
    for (unsigned int vertex = 0; vertex < firstDirectedEdge.size(); vertex++)
        outFile << "FirstDirectedEdge " << vertex<< "   " << std::fixed << firstDirectedEdge[vertex] << std::endl;

    // write and the faces - increment is taken care of internally
	int faceId = 0;
    for (unsigned int face = 0; face < faceVertices.size(); face+=3)
        { // per face
        outFile << "Face  " << faceId++ << "   ";
        
        // write three vertices
        outFile << faceVertices[face] << "  ";
        outFile << faceVertices[face+1] << "  ";
        outFile << faceVertices[face+2];
            
        outFile << std::endl;
        } // per face

	// write the other halves
	for (unsigned int dirEdge = 0; dirEdge < otherHalf.size(); dirEdge++)
		outFile << "OtherHalf " << dirEdge << "   " << otherHalf[dirEdge] << std::endl;

	outFile.close();
}