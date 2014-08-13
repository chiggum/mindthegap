//bitmap.cpp
#include "bitmap.h"
#include "fstream"
#include <iostream>
#include "debug.h"
#include <string>
#include <algorithm>
#include <cmath>
#include "medianBlur.h"
#include "posterize.h"
#include "info.h"

//Constructor
Bitmap::Bitmap(class Info *info)
{
	inputParam = info;
	uint h, w;
	std::vector<uchar> v;	//the raw pixel

	decodeOneStep(inputParam->inputFileName.c_str(), h, w, v);
	//the pixels are now in the vector "v", 4 bytes per pixel, ordered RGBARGBA

	//Initialization and allocation of private variables
	intialize(&orig, h, w);
	intialize(&pre, h + 2, w + 2);
	intialize(&pop, 2 * (h + 2), 2 * (w + 2));

	pixToNodeMap = new int*[2 * (h + 2)];
	for(uint i = 0; i < 2 * (h + 2); ++i)
	{
		pixToNodeMap[i] = new int[2 * (w + 2)];
	}

	numControlPoints = 0;
	codedImage = NULL;
	graph = NULL;

	//Assigning values to the pixels of original image matrix
	for(uint i = 0; i < h; ++i)
	{
		for(uint j = 0; j < w; ++j)
		{
			orig->pixMap[i][j].r = v[4 * (i * w + j)];
			orig->pixMap[i][j].g = v[4 * (i * w + j) + 1];
			orig->pixMap[i][j].b = v[4 * (i * w + j) + 2];
		}
	}

	//borderPixel initialization
	if(inputParam->bgColorProvided)
	{
		borderPixel = inputParam->bgColor;
	}
	else
	{
		//default value of border pixel if bgColor not provided
		/*
		 * let F(pixel) = pixel.r + 256 * pixel.g + 256 * 256 * pixel.b
		 * step 1: calculate F(.) for four corner pixels of original image and store in an array
		 * step 2: sort the array in ascending order
		 * step 3: find median of four values i.e. mean of 2nd and 3rd value
		 * step 4: interpolate the median to get r, g, b values.
		 */
		uint a[4];
		a[0] = orig->pixMap[0][0].r 		+ orig->pixMap[0][0].g * 256 		 + orig->pixMap[0][0].b * 256 * 256;
		a[1] = orig->pixMap[0][w - 1].r 	+ orig->pixMap[0][w - 1].g * 256 	 + orig->pixMap[0][w - 1].b * 256 * 256;
		a[2] = orig->pixMap[h - 1][0].r 	+ orig->pixMap[h - 1][0].g * 256 	 + orig->pixMap[h - 1][0].b * 256 * 256;
		a[3] = orig->pixMap[h - 1][w - 1].r + orig->pixMap[h - 1][w - 1].g * 256 + orig->pixMap[h - 1][w - 1].b * 256 * 256;

		std::sort(a, a + 4);
		double med = (a[1] + a[2]) * 0.5;
		uint median = floor(med);

		//interpolation
		borderPixel.b = median / (256 * 256);
		median = median - borderPixel.b * 256 * 256;
		borderPixel.g = median / 256;
		median = median - borderPixel.g * 256;
		borderPixel.r = median;
	}

/******************************************************************/
/***********************TEST 1:DEBUGGING MODE**********************/
/******************************************************************/

#ifdef _TEST_1_ 

	//file to store r g b values of each pixel in row major order under test 1
	std::string path = ROOT_DIR;
	std::ofstream ofsTest1((path + "/check/test_1_cpp.txt").c_str(), std::ofstream::out);
	ofsTest1 << h << std::endl;
	ofsTest1 << w << std::endl;
	
	for(uint i = 0; i < h; ++i)
	{
		for(uint j = 0; j < w; ++j)
		{
			//print r g b int values in new lines
			ofsTest1 << (uint)orig->pixMap[i][j].r << std::endl; 
			ofsTest1 << (uint)orig->pixMap[i][j].g << std::endl;
			ofsTest1 << (uint)orig->pixMap[i][j].b << std::endl;
		}
	}

	//close test_1_cpp.txt file
	ofsTest1.close();

#endif	

/******************************************************************/
/***********************TEST 1:DEBUGGING MODE END******************/
/******************************************************************/
}

//Decode from disk to raw pixels with a single function call
void Bitmap::decodeOneStep(const char* filename, uint &height, uint &width, std::vector<uchar> &image)
{
  //decode
  unsigned error = lodepng::decode(image, width, height, filename);

  //if there's an error, display it
  if(error) 
  	std::cout << "decoder error " << error << ": " << lodepng_error_text(error) << std::endl;
}

//Initialize ImageMatrix
void Bitmap::intialize(ImageMatrix **m, uint h, uint w)
{
	*m = new ImageMatrix;
	(*m)->height = h;
	(*m)->width = w;
	(*m)->pixMap = new pixel*[h];
	for(uint i = 0; i < h; ++i)
	{
		(*m)->pixMap[i] = new pixel[w];
	}
}

//Destructor
Bitmap::~Bitmap()
{
	//delete pixToNodeMap
	for(uint i = 0; i < pop->height; ++i)
	{
		delete[] pixToNodeMap[i];
	}
	delete[] pixToNodeMap;

	//delete orig
	del(orig);

	//delete pre
	del(pre);

	//delete pop
	del(pop);

	//delete coded image
	delete codedImage;

	//delete graph
	delete graph;
}

//Delete an ImageMatrix
void Bitmap::del(ImageMatrix *a)
{
	for(uint i = 0; i < a->height; ++i)
	{
		delete[] a->pixMap[i];
	}
	delete[] a->pixMap;
	delete a;
}

/*
 * Processes the orig image
 */
void Bitmap::processImage()
{
	//if noisy switch is ON
	if(inputParam->switchNoisy)
		removeNoise();

	preprocess();
	popoutBoundaries();

	//create codeImage object corresponding to popped out bounadires bitmap
	codedImage = new CodeImage(pop);

	//create graph object corresponding to popped out boundaries image
	graph = new Graph(pop->height, pop->width, *codedImage->getColCode(), pop);

	//detect control points and initialize graph with boundary pixels as vertices
	detectControlPoints();

	//form adjacency list of graph
	formAdjacencyList();

	//removes the connections defined Dangerous(see graph::formLineSegments)
	removeDangerousConnections();

	//forms line segments from control point to control point
	graph->formLineSegments(1);

	//forms island line segments
	graph->formLineSegments(2);

	//assign the curves id/num to regions of pop bitmap
	graph->assignCurveNumToRegion();

	//preprocess the line segments before passing to fitting routine
	graph->preprocessLineSegments();

	//form curves by passing series of points on line segments to fitting routine
	graph->formCurves(inputParam->toleranceCurve, inputParam->toleranceLine);

	//form boundaries/paths surrounding each region
	graph->processRegions();
}

/*
 * Remove noise from image
 * step 1: perform median blur with user provided kernel size(default 3X3).
 * step 2: perform posterization on blurred image with user provided number
 * 		   of clusters and max iterations (default= 8 and 50 respectively).
 * step 3: remove turds of size <= turdSize provided (default=10). This will
 * 		   reduce no. of control points
 */
void Bitmap::removeNoise()
{
	medianBlur(inputParam->medianBlurKernelSize, orig);

	posterize(orig, inputParam->numClusters, inputParam->KMeansMaxIters);

	//create a coded image corresponding to orig i.e. posterized image (here)
	CodeImage *myCodedImage = new CodeImage(orig);

	//process regions i.e. remove turds of size <= turdSize
	myCodedImage->processRegions(inputParam->turdSize);

	ImageMatrix *m = new ImageMatrix[1];
	m = myCodedImage->getFinalImage();

	//copy back final non-noisy image to orig
	for(uint i = 0; i < orig->height; ++i)
	{
		for(uint j = 0; j < orig->width; ++j)
		{
			orig->pixMap[i][j] = m->pixMap[i][j];
		}
	}

	delete myCodedImage;
	delete m;
}

//Write output svg in outFileName by moving control to writeOutput method of graph
void Bitmap::writeOuputSVG()
{
	graph->writeOuput(inputParam->outFileName, borderPixel);
}

/*
 * Appends row at top and bottom and col at left and right of 
 * orig image with borderPixel
 * 
 * Also, finds a unique rgb color to be used as popped out boundary color.
 * The method to find unique color assumes that there exit atleast a single
 * value in [0, 255] which is not taken by any pixel's r or g or b channel.
 * For example: if no pixel of orig bitmap has r channel == 100 then a unique
 * pixel would be (100, 255, 255) which is then assigned to boundary pixel.
 */
void Bitmap::preprocess()
{
	uint h = pre->height;
	uint w = pre->width;
	int M[3][256] = {{0}};

	for(uint i = 0; i < h; ++i)
	{
		for(uint j = 0; j < w; ++j)
		{
			if(i == 0 || j == 0 || i == h - 1 || j == w - 1)
			{
				pre->pixMap[i][j] = borderPixel;
			}
			else
			{
				pre->pixMap[i][j] = orig->pixMap[i - 1][j - 1];
				M[0][orig->pixMap[i - 1][j - 1].r] = M[0][orig->pixMap[i - 1][j - 1].g] = M[0][orig->pixMap[i - 1][j - 1].b] = -1;				
			}
		}
	}

	//Getting unique boundaryPixel
	bool outOfLoop = false;
	for(uint i = 0; i < 3; ++i)
	{
		for(uint j = 0; j < 256; ++j)
		{
			if(M[i][j] != -1)
			{
				switch(i)
				{
					case 0:
						boundaryPixel.r = j;
						boundaryPixel.g = boundaryPixel.b = 255;
						break;
					case 1:
						boundaryPixel.g = j;
						boundaryPixel.r = boundaryPixel.b = 255;
						break;
					case 2:
						boundaryPixel.b = j;
						boundaryPixel.g = boundaryPixel.r = 255;
						break;
				}
				outOfLoop = true;
				break;
			}
		}
		if(outOfLoop)
			break;
	}
}

/* 
 * Popping Out Boundaries.
 * 	Suppose the input image is:
 * 		orig =
 * 			1	1	2
 * 			1	2	2
 * 			3 	3	3
 * 	where 1, 2 and 3 represents colors/pixels.
 * 	Let the background color be 0.
 * 
 * 	Then the preprocessed image will be:
 * 		pre =
 * 			0	0	0	0	0
 * 			0	1	1	2	0
 * 			0	1	2	2	0
 * 			0	3	3	3	0
 * 			0	0	0	0	0
 * 
 * 	popoutBoundaries() routine:
 * 
 * 		Let pop be our final popped out boundaries image.
 * 		step 1: initializes pop with boundary pixels. After intialization,
 * 			pop = 
 * 				0	0	0	0	0	0	0	0	0	0
 * 				0	0	0	0	0	0	0	0	0	0
 * 				0	0	0	0	0	0	0	0	0	0
 * 				0	0	0	0	0	0	0	0	0	0
 * 				0	0	0	0	0	0	0	0	0	0
 * 				0	0	0	0	0	0	0	0	0	0
 * 				0	0	0	0	0	0	0	0	0	0
 * 				0	0	0	0	0	0	0	0	0	0
 * 				0	0	0	0	0	0	0	0	0	0
 * 				0	0	0	0	0	0	0	0	0	0
 * 
 * 		step 2: Assuming ZERO INDEXED rows and columns
 * 			for i goes from 1st row to last row of pre
 * 				focus_row := (2 * i + 1)th row
 * 				for j goes from 1st col to last col of pre
 * 					focus_col := (2 * j + 1)th col
 * 
 * 					pixel_in_pop_at(focus_row, focus_col) := pixel_in_pre_at(i,j) 
 * 
 * 					if(pixel_in_pre_at(i,j) == pixel_in_pre_at(i, j-1))
 * 						pixel_in_pop_at(focus_row, focus_col-1) := pixel_in_pre_at(i,j-1)
 * 					else
 * 						pixel_in_pop_at(focus_row, focus_col-1) := boundary pixel
 * 					endif
 * 
 * 					if(pixel_in_pre_at(i,j) == pixel_in_pre_at(i-1, j-1))
 * 						pixel_in_pop_at(focus_row-1, focus_col-1) := pixel_in_pre_at(i-1,j-1)
 * 					else
 * 						pixel_in_pop_at(focus_row-1, focus_col-1) := boundary pixel
 * 					endif
 * 
 * 					if(pixel_in_pre_at(i,j) == pixel_in_pre_at(i-1, j))
 * 						pixel_in_pop_at(focus_row-1, focus_col) := pixel_in_pre_at(i-1,j)
 * 					else
 * 						pixel_in_pop_at(focus_row-1, focus_col) := boundary pixel
 * 					endif
 * 				endfor
 * 			endfor
 * 
 * 			Using the above algorithm, (X represents boundary pixels)
 * 				pop =
 * 					0	0	0	0	0	0	0	0	0	0
 * 					0	0	0	0	0	0	0	0	0	0
 * 					0	0	X	X	X	X	X	X	0	0
 * 					0	0	X	1	1	1	X	2	X	0
 * 					0	0	X	1	X	X	X	2	X	0
 * 					0	0	X	1	X	2	2	2	X	0
 * 					0	0	X	X	X	X	X	X	X	0
 * 					0	0	X	3	3	3	3	3	X	0
 * 					0	0	0	X	X	X	X	X	X	0
 * 					0	0	0	0	0	0	0	0	0	0
 */
void Bitmap::popoutBoundaries()
{
	uint h = pop->height;
	uint w = pop->width;
	
	//step 1: initialization
	for(uint i = 0; i < h; ++i)
	{
		for(uint j = 0; j < w; ++j)
		{
			pop->pixMap[i][j] = borderPixel;
		}
	}

	//step 2: popping out boundaries
	for(uint i = 1; i < h * 0.5; ++i)
	{
		uint focus_i = 2 * i + 1;
		for(uint j = 1; j < w * 0.5; ++j) 
		{
			uint focus_j = 2 * j + 1;

			pop->pixMap[focus_i][focus_j] = pre->pixMap[i][j];

			if(!ifEqualPixel(pre->pixMap[i][j], pre->pixMap[i][j - 1]))
			{
				pop->pixMap[focus_i][focus_j - 1] = boundaryPixel;
			}
			else
			{
				pop->pixMap[focus_i][focus_j - 1] = pre->pixMap[i][j - 1];
			}
			if(!ifEqualPixel(pre->pixMap[i][j], pre->pixMap[i - 1][j - 1]))
			{
				pop->pixMap[focus_i - 1][focus_j - 1] = boundaryPixel;
			}
			else
			{
				pop->pixMap[focus_i - 1][focus_j - 1] = pre->pixMap[i - 1][j - 1];
			}
			if(!ifEqualPixel(pre->pixMap[i][j], pre->pixMap[i - 1][j]))
			{
				pop->pixMap[focus_i - 1][focus_j] = boundaryPixel;
			}
			else
			{
				pop->pixMap[focus_i - 1][focus_j] = pre->pixMap[i - 1][j];
			}
		}
	}

/******************************************************************/
/***********************TEST 2:DEBUGGING MODE**********************/
/******************************************************************/

#ifdef _TEST_2_

	std::string path = ROOT_DIR;
	//file to store r g b values of each pixel in row major order under test 2
	std::ofstream ofsTest2((path + "/check/test_2_cpp.txt").c_str(), std::ofstream::out);
	ofsTest2 << h << std::endl;
	ofsTest2 << w << std::endl;

	for(uint i = 0; i < h; ++i)
	{
		for(uint j = 0; j < w; ++j)
		{
			//print r g b int values in new lines
			ofsTest2 << (int)pop->pixMap[i][j].r << std::endl;
			ofsTest2 << (int)pop->pixMap[i][j].g << std::endl;
			ofsTest2 << (int)pop->pixMap[i][j].b << std::endl;
		}
	}

	//close test_2_cpp.txt file
	ofsTest2.close();

#endif

/******************************************************************/
/***********************TEST 2:DEBUGGING MODE END******************/
/******************************************************************/

/******************************************************************/
/***********************EVAL 1:DEBUGGING MODE**********************/
/******************************************************************/

#ifdef _EVAL_1_

	//generate some image
	#ifndef _TEST_2_
		//prevent double declaration
		std::string path = ROOT_DIR;
	#endif
	std::vector<unsigned char> image;
	image.resize(pop->width * pop->height * 4);
	for(uint y = 0; y < pop->height; y++)
	{
		for(uint x = 0; x < pop->width; x++)
		{
			image[4 * pop->width * y + 4 * x + 0] = pop->pixMap[y][x].r;
			image[4 * pop->width * y + 4 * x + 1] = pop->pixMap[y][x].g;
			image[4 * pop->width * y + 4 * x + 2] = pop->pixMap[y][x].b;
			image[4 * pop->width * y + 4 * x + 3] = 255;
		}
	}

	encodeOneStep((path + "/check/eval_1_popout.png").c_str(), image, pop->width, pop->height);
	image.clear();

#endif

/******************************************************************/
/***********************EVAL 1:DEBUGGING MODE END******************/
/******************************************************************/

}

/*
 * Detects control points which are boundary pixels where 3 or more regions meet.
 * For example: the middle X is a control point in (1) not in (2).
 *	1 X 2 	1 X 2
 *  X X 2 	1 X 2
 *  3 3 X 	1 X 2
 *   (1)     (2)
 * Also, initializes the graph with its vertices as boundary pixels.
*/
void Bitmap::detectControlPoints()
{
	uint h = pop->height;
	uint w = pop->width;
	uint vertexCounter = 0;
	std::vector<pixel> temp;

	for(uint i = 1; i < h - 1; ++i)
	{
		for(uint j = 1; j < w - 1; ++j)
		{
			//if pixel at (i,j) is a boundary pixel
			if(ifEqualPixel(pop->pixMap[i][j], boundaryPixel))
			{
				//assign vertexCounter value to (i,j)th entry of pixToNodeMap
				pixToNodeMap[i][j] = vertexCounter;

				//add vertex in graph corresponding to this boundary pixel(i,j)
				graph->addVertex(vertexCounter, i, j, 0, 0, 0);

				temp.clear();

				//if pixel at (i,j-1) is a unique (newly found) adjacent region pixel
				if(checkUniqueRegionPixel(pop->pixMap[i][j - 1], temp))
					//add adjacent region code to adjacent region vector of vertex in graph
					graph->addRegion(vertexCounter, codedImage->getMatrix()[i][j - 1]);	

				if(checkUniqueRegionPixel(pop->pixMap[i][j + 1], temp))
					graph->addRegion(vertexCounter, codedImage->getMatrix()[i][j + 1]);

				if(checkUniqueRegionPixel(pop->pixMap[i - 1][j - 1], temp))
					graph->addRegion(vertexCounter, codedImage->getMatrix()[i - 1][j - 1]);

				if(checkUniqueRegionPixel(pop->pixMap[i - 1][j + 1], temp))
					graph->addRegion(vertexCounter, codedImage->getMatrix()[i - 1][j + 1]);

				if(checkUniqueRegionPixel(pop->pixMap[i + 1][j - 1], temp))
					graph->addRegion(vertexCounter, codedImage->getMatrix()[i + 1][j - 1]);

				if(checkUniqueRegionPixel(pop->pixMap[i + 1][j + 1], temp))
					graph->addRegion(vertexCounter, codedImage->getMatrix()[i + 1][j + 1]);

				if(checkUniqueRegionPixel(pop->pixMap[i - 1][j], temp))
					graph->addRegion(vertexCounter, codedImage->getMatrix()[i - 1][j]);

				if(checkUniqueRegionPixel(pop->pixMap[i + 1][j], temp))
					graph->addRegion(vertexCounter, codedImage->getMatrix()[i + 1][j]);

				//if more than or equal to 3 different adjacent pixel
				if(temp.size() >= 3)
				{
					//set node with index = vertexCounter, a control point in graph
					graph->setControlPoint(vertexCounter, 1);
				}

				vertexCounter++;
			}
			else
			{
				//not a white pixel
				pixToNodeMap[i][j] = -1;
				if(ifEqualPixel(pop->pixMap[i][j-1],boundaryPixel) && ifEqualPixel(pop->pixMap[i][j+1], boundaryPixel) && ifEqualPixel(pop->pixMap[i-1][j], boundaryPixel) && ifEqualPixel(pop->pixMap[i+1][j], boundaryPixel))
				{
					pixToNodeMap[i][j] = -2;	//danger pixel
				}
			}
		}
	}

	temp.clear();

/******************************************************************/
/***********************EVAL 2:DEBUGGING MODE**********************/
/******************************************************************/

#ifdef _EVAL_2_
	
	std::string path = ROOT_DIR;
	std::vector<unsigned char> image;
	image.resize(w * h * 4);

	for(uint i = 1; i < h - 1; ++i)
	{
		for(uint j = 1; j < w - 1; ++j)
		{

			image[4 * w * i + 4 * j + 0] = pop->pixMap[i][j].r;
			image[4 * w * i + 4 * j + 1] = pop->pixMap[i][j].g;
			image[4 * w * i + 4 * j + 2] = pop->pixMap[i][j].b;
			image[4 * w * i + 4 * j + 3] = 255;

			if(pixToNodeMap[i][j] != -1 && graph->checkIfCntrlPt(pixToNodeMap[i][j]))
			{
				image[4 * w * i + 4 * j + 0] = 0;
				image[4 * w * i + 4 * j + 1] = 0;
				image[4 * w * i + 4 * j + 2] = 255;
				image[4 * w * i + 4 * j + 3] = 255;
			}
			else if(pixToNodeMap[i][j] == -2)
			{
				image[4 * w * i + 4 * j + 0] = 255;
				image[4 * w * i + 4 * j + 1] = 0;
				image[4 * w * i + 4 * j + 2] = 0;
				image[4 * w * i + 4 * j + 3] = 255;
			}
		}
	}

	encodeOneStep((path + "/check/eval_2_controlPt.png").c_str(), image, w, h);
	image.clear();

#endif

/******************************************************************/
/***********************EVAL 2:DEBUGGING MODE END******************/
/******************************************************************/

/******************************************************************/
/***********************TEST 4:DEBUGGING MODE**********************/
/******************************************************************/

#ifdef _TEST_4_
	
	#ifndef _EVAL_2_
		std::string path = ROOT_DIR;
	#endif
	//file to store coordinate values of control points under test 4
	std::ofstream ofsTest4((path + "/check/test_4_cpp.txt").c_str(), std::ofstream::out);
	ofsTest4 << h << std::endl;
	ofsTest4 << w << std::endl;
	for(uint i = 1; i < h - 1; ++i)
	{
		for(uint j = 1; j < w - 1; ++j)
		{
			if(pixToNodeMap[i][j] != -1 && graph->checkIfCntrlPt(pixToNodeMap[i][j]))
			{
				//print coordinates of control point
				ofsTest4 << i << std::endl;
				ofsTest4 << j << std::endl;
			}
		}
	}

	//close file
	ofsTest4.close();

#endif

/******************************************************************/
/***********************TEST 4:DEBUGGING MODE END******************/
/******************************************************************/

}

/*
 * if(pixel "a" is not a boundary pixel and is not contained in vector v)
 *     then push a in v and return 1 
 * else return 0.
 * vector v contains all the region pixels adjacent to a node(= boundarypixel),
 * detected before calling this routine.
 */
int Bitmap::checkUniqueRegionPixel(pixel a, std::vector<pixel> &v)
{
	int ret = 0;
	//if "a" is not a boundary pixel
	if(!ifEqualPixel(a, boundaryPixel))
	{
		uint flag = 0;
		for(uint i = 0; i < v.size(); ++i)
		{
			//if "a" found in v
			if(ifEqualPixel(a, v[i]))
			{
				flag = 1;
				break;
			}
		}

		//if "a" not found in v
		if(flag == 0)
		{
			v.push_back(a);
			ret = 1;
		}
	}

	return ret;
}

/*
 * Forms adjacency list in graph.
 * Connects two vertices in the graph if they are adjacent.
 * The connection is bidirectional(undirected garph).
 * Since only boundary pixels are nodes/vertices, we loop over all the pixels and
 * check if the pixel is a boundary pixel,
 * if yes: then connects this pixel/vertex with any adjacent boundary pixel/vertex.
 */
void Bitmap::formAdjacencyList()
{
	uint h = pop->height;
	uint w = pop->width;

	for(uint i = 1; i < h - 1; ++i)
	{
		for(uint j = 1; j < w - 1; ++j)
		{
			//if pixel at (i,j) is a boundary pixel
			if(ifEqualPixel(pop->pixMap[i][j], boundaryPixel))
			{
				//if pixel at (i-1,j) is a boundary pixel
				if(ifEqualPixel(pop->pixMap[i - 1][j], boundaryPixel))
					//connect vertex at (i,j) and at (i,j-1)
					graph->addEdge(pixToNodeMap[i][j], pixToNodeMap[i - 1][j]);

				if(ifEqualPixel(pop->pixMap[i + 1][j], boundaryPixel))
					graph->addEdge(pixToNodeMap[i][j], pixToNodeMap[i + 1][j]);

				if(ifEqualPixel(pop->pixMap[i - 1][j - 1], boundaryPixel))
					graph->addEdge(pixToNodeMap[i][j], pixToNodeMap[i - 1][j - 1]);

				if(ifEqualPixel(pop->pixMap[i - 1][j + 1], boundaryPixel))
					graph->addEdge(pixToNodeMap[i][j], pixToNodeMap[i - 1][j + 1]);

				if(ifEqualPixel(pop->pixMap[i + 1][j - 1], boundaryPixel))
					graph->addEdge(pixToNodeMap[i][j], pixToNodeMap[i + 1][j - 1]);

				if(ifEqualPixel(pop->pixMap[i + 1][j + 1], boundaryPixel))
					graph->addEdge(pixToNodeMap[i][j], pixToNodeMap[i + 1][j + 1]);

				if(ifEqualPixel(pop->pixMap[i][j - 1], boundaryPixel))
					graph->addEdge(pixToNodeMap[i][j], pixToNodeMap[i][j - 1]);

				if(ifEqualPixel(pop->pixMap[i][j + 1], boundaryPixel))
					graph->addEdge(pixToNodeMap[i][j], pixToNodeMap[i][j + 1]);
			}
		}
	}
}

/*
 * Removes dangerous connections around dangerous points.
 * Please refer graph::formLineSegments to know about dangerous points and connections
 */
void Bitmap::removeDangerousConnections()
{
	uint h = pop->height;
	uint w = pop->width;

	for(uint i = 1; i < h - 1; ++i)
	{
		for(uint j = 1; j < w - 1; ++j)
		{
			if(pixToNodeMap[i][j] == -2)
			{
				if(!graph->checkIfUsedUp(pixToNodeMap[i][j-1]) && !graph->checkIfUsedUp(pixToNodeMap[i][j+1]) && !graph->checkIfUsedUp(pixToNodeMap[i-1][j]) && !graph->checkIfUsedUp(pixToNodeMap[i+1][j]))
				{
					pixel midPixel = pop->pixMap[i][j];
					if(ifEqualPixel(pop->pixMap[i-1][j-1], midPixel) && ifEqualPixel(pop->pixMap[i+1][j+1], midPixel))
					{
						graph->removeConnection(pixToNodeMap[i-1][j], pixToNodeMap[i][j-1]);
						graph->removeConnection(pixToNodeMap[i+1][j], pixToNodeMap[i][j+1]);
					}
					else if(ifEqualPixel(pop->pixMap[i+1][j-1], midPixel) && ifEqualPixel(pop->pixMap[i-1][j+1], midPixel))
					{
						graph->removeConnection(pixToNodeMap[i-1][j], pixToNodeMap[i][j+1]);
						graph->removeConnection(pixToNodeMap[i+1][j], pixToNodeMap[i][j-1]);
					}
				}
			}
		}
	}
}