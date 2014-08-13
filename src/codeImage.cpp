//codeImage.cpp
#include "codeImage.h"
#include "debug.h"
#include <fstream>
#include <string>
#include "bitmap.h"
#include <queue>
#include <utility>
#include "GraphicsGems.h"
#include <cstdlib>
#include <vector>
#include <iostream>
#include "util.h"

//Constructor
CodeImage::CodeImage(ImageMatrix* m)
{
	height = m->height;
	width = m->width;

	//initialize mat
	mat = new int*[height];
	for(uint i = 0; i < height; ++i)
	{
		mat[i] = new int[width];
		for(uint j = 0; j < width; ++j)
			mat[i][j] = -1;	//-1 represents not coded yet
	}

	codeImage(m);
	numDisjointRegions = colorCode.size();


/******************************************************************/
/***********************TEST 3:DEBUGGING MODE**********************/
/******************************************************************/

#ifdef _TEST_3_

	//generate some image
	std::string path = ROOT_DIR;
	std::vector<unsigned char> image;
	image.resize(width * height * 4);
	for(uint y = 0; y < height; y++)
	{
		for(uint x = 0; x < width; x++)
		{
			image[4 * width * y + 4 * x + 0] = colorCode[mat[y][x]].r;
			image[4 * width * y + 4 * x + 1] = colorCode[mat[y][x]].g;
			image[4 * width * y + 4 * x + 2] = colorCode[mat[y][x]].b;
			image[4 * width * y + 4 * x + 3] = 255;
		}
	}

	encodeOneStep((path + "/check/test_3_revertCoded.png").c_str(), image, width, height);
	image.clear();

#endif

/******************************************************************/
/***********************TEST 3:DEBUGGING MODE END******************/
/******************************************************************/
}

//Destructor
CodeImage::~CodeImage()
{
	//delete mat
	for(uint i = 0; i < height; ++i)
	{
		delete[] mat[i];
	}
	delete[] mat;

	//clear colorCode
	colorCode.clear();
}

/*
 * codeImage.
 * Codes the input image using operation similar to flood fill
 *
 * 	Let m be the input ImageMatrix:
 * 		m = 
 * 			*	*	*	*	*	*	*	*	*	*
 *			*	*	*	*	*	*	*	*	*	*
 *			*	*	#	#	#	#	#	#	*	*
 *			*	*	#	@	@	@	#	$	#	*
 *			*	*	#	@	#	#	#	$	#	*
 *			*	*	#	@	#	$	$	$	#	*
 *			*	*	#	#	#	#	#	#	#	*
 *			*	*	#	%	%	%	%	%	#	*
 *			*	*	*	#	#	#	#	#	#	*
 *			*	*	*	*	*	*	*	*	*	*
 *
 *	Let mat be the matrix that stores the code corresponding to each pixel.
 *	Note that: intially mat[i][j] = -1 for all i and j.
 *
 *	A high level overview of algorithm: Assuming ZERO INDEXED rows and cols
 *		counter = 0;
 *		for i = 0 to last row
 *			for j = 0 to last row
 *				if((i,j) is not yet coded)
 *					code pixel(i,j) with counter i.e. mat[i][j] = counter 
 *					Q be an empty queue.
 *					push pixel(i, j) in Q.
 *					while(Q is not empty)
 *						x = Q.pop().
 *						assign all pixels adjacent to x that are equal to x, the same code i.e. counter and push them in queue Q
 *					endwhile
 *					increment counter by 1.
 *				endif
 *			endfor
 *		endfor
 *
 *	For the particular m that we chose, the number of times the routine enters into if statement is 5.
 *	After it comes out 1st time:
 *		0	0	0	0	0	0	0	0	0	0
 *		0	0	0	0	0	0	0	0	0	0
 *		0	0	#	#	#	#	#	#	0	0
 *		0	0	#	@	@	@	#	$	#	0
 *		0	0	#	@	#	#	#	$	#	0
 *		0	0	#	@	#	$	$	$	#	0
 *		0	0	#	#	#	#	#	#	#	0
 *		0	0	#	%	%	%	%	%	#	0
 *		0	0	0	#	#	#	#	#	#	0
 *		0	0	0	0	0	0	0	0	0	0
 *
 *	After it comes out 2nd time:
 *		0	0	0	0	0	0	0	0	0	0
 *		0	0	0	0	0	0	0	0	0	0
 *		0	0	1	1	1	1	1	1	0	0
 *		0	0	1	@	@	@	1	$	1	0
 *		0	0	1	@	1	1	1	$	1	0
 *		0	0	1	@	1	$	$	$	1	0
 *		0	0	1	1	1	1	1	1	1	0
 *		0	0	1	%	%	%	%	%	1	0
 *		0	0	0	1	1	1	1	1	1	0
 *		0	0	0	0	0	0	0	0	0	0
 *
 *	After it comes out 3rd time:
 *		0	0	0	0	0	0	0	0	0	0
 *		0	0	0	0	0	0	0	0	0	0
 *		0	0	1	1	1	1	1	1	0	0
 *		0	0	1	2	2	2	1	$	1	0
 *		0	0	1	2	1	1	1	$	1	0
 *		0	0	1	2	1	$	$	$	1	0
 *		0	0	1	1	1	1	1	1	1	0
 *		0	0	1	%	%	%	%	%	1	0
 *		0	0	0	1	1	1	1	1	1	0
 *		0	0	0	0	0	0	0	0	0	0
 *
 *	After it comes out 4th time:
 *		0	0	0	0	0	0	0	0	0	0
 *		0	0	0	0	0	0	0	0	0	0
 *		0	0	1	1	1	1	1	1	0	0
 *		0	0	1	2	2	2	1	3	1	0
 *		0	0	1	2	1	1	1	3	1	0
 *		0	0	1	2	1	3	3	3	1	0
 *		0	0	1	1	1	1	1	1	1	0
 *		0	0	1	%	%	%	%	%	1	0
 *		0	0	0	1	1	1	1	1	1	0
 *		0	0	0	0	0	0	0	0	0	0
 *
 *	After it comes out 5th time:
 *		0	0	0	0	0	0	0	0	0	0
 *		0	0	0	0	0	0	0	0	0	0
 *		0	0	1	1	1	1	1	1	0	0
 *		0	0	1	2	2	2	1	3	1	0
 *		0	0	1	2	1	1	1	3	1	0
 *		0	0	1	2	1	3	3	3	1	0
 *		0	0	1	1	1	1	1	1	1	0
 *		0	0	1	3	3	3	3	3	1	0
 *		0	0	0	1	1	1	1	1	1	0
 *		0	0	0	0	0	0	0	0	0	0
 *
 * Also, regionPixelCoord[i] stores all the pixels' coordinates
 * corresponding to ith region.
 */
void CodeImage::codeImage(ImageMatrix *m)
{
	uint counter = 0;
	for(uint ix = 0; ix < height; ++ix)
	{
		for(uint jy = 0; jy < width; ++jy)
		{
			//if pixel at (ix,jy) is not yet coded
			if(mat[ix][jy] == -1)
			{
				//a temporary vector that stores the coordinates
				//of all the pixels in the region with code counter
				std::vector<Point2> temp;
				Point2 *tempCoord = new Point2[1];
				tempCoord->x = ix;
				tempCoord->y = jy;
				temp.push_back(*tempCoord);

				mat[ix][jy] = counter;

				//storing the color of this region in colorCode
				pixel *myPixel = new pixel[1];
				myPixel->r = m->pixMap[ix][jy].r;
				myPixel->g = m->pixMap[ix][jy].g;
				myPixel->b = m->pixMap[ix][jy].b;
				colorCode.push_back(*myPixel);

				//an empty queue which hold integer coordinates of pixel
				std::queue< std::pair<int, int> >Q;
				
				//push (ix,jy) in Q
				Q.push(std::make_pair(ix,jy));
				while(!Q.empty())
				{
					int i = Q.front().first;
					int j = Q.front().second;

					Q.pop();

					//checking all eight coordinates surrounding (i,j) with boundary conditions checked
					for(int l = i - 1; l <= i + 1; ++l)
					{
						for(int t = j - 1; t <= j + 1; ++t)
						{
							if(l >= 0 && l <= (int)height - 1 && t >= 0 && t <= (int)width - 1 && (l != i || t != j))
							{
								//if pixel at (l,t) is not coded yet and pixel at (l,t) == (i,j)
								if(mat[l][t] == -1 && ifEqualPixel(m->pixMap[i][j], m->pixMap[l][t]))
								{
									mat[l][t] = counter;
									tempCoord = new Point2[1];
									tempCoord->x = l;
									tempCoord->y = t;

									//push the coordinate in temp and Q
									temp.push_back(*tempCoord);
									Q.push(std::make_pair(l,t));
								}
							}
						}
					}
				}
				//push all the coordinates of the pixels in the region with code counter in regionPixelCoord
				regionPixelCoord.push_back(temp);
				temp.clear();
				counter = counter + 1;
			}
		}
	}
}

//Returns the coded image matrix
int** CodeImage::getMatrix()
{
	return mat;
}

//Returns the vector of colors of unique regions
std::vector<pixel>* CodeImage::getColCode()
{
	return &colorCode;
}

/*
 * Processes regions to remove regions with size <= threshold
 * Here, size of a region = no. of pixels in that region
 */
void CodeImage::processRegions(uint threshold)
{
	//five passes are sufficient to remove all the major turds
	int pass = 5;
	while(pass > 0)
	{
		for(uint i = 0; i < regionPixelCoord.size(); ++i)
		{
			//if no. of pixel in the region <= threshold
			if(regionPixelCoord[i].size() <= threshold)
			{
				//dissolve/merge region with surrounding regions
				dissolveRegion(i);
			}
		}
		pass--;
	}
}

/*
 *  dissolveRegion.
 *  It dissolves a region in its surrounding regions.
 *  A high level overview of the algorithm:
 *  	Suppose rgn be the region that needs to be dissolved/merged.
 *  		while(rgn contains a pixel i.e. it is not empty)
 *  			iterate over all the pixels of the region
 *  				let boundElem and boundCount be two empty vectors.
 *					iterate over all surrounding pixels of this pixel
 *  					if(this surrounding pixel has a different code than this pixel)
 *  						if(this different code pixel found is already in boundElem)
 *  							increase boundCount corresponding to this different pixel
 *  						else
 *  							push this different pixel in boundElem 
 *								and initialize boundCount corresponding to this different pixel with 1
 *  						endif
 *  					endif
 *					enditerate
 *  			enditerate
 *  			if(boundElem is not empty)
 *  				//This means that this pixel is a boundary pixel
 *  				find the boundElem which has highest boundCount and assign it's code to this pixel.
 *  				push this pixel's coordinate to the region which is assigned to this pixel.
 *					Also, erase this pixel from rgn.
 *  			endif
 *  		endwhile
 *
 */
void CodeImage::dissolveRegion(int ind)
{
	while(!regionPixelCoord[ind].empty())
	{
		//iterating over all pixels of the region
		for(uint i = 0; i < regionPixelCoord[ind].size(); ++i)
		{
			int ix = regionPixelCoord[ind][i].x;
			int jy = regionPixelCoord[ind][i].y;

			std::vector<int> boundElem, boundCount;

			//iterate over surrounding pixel
			for(int l = ix - 1; l <= ix + 1; ++l)
			{
				for(int m = jy - 1; m <= jy + 1; ++m)
				{
					if(l >= 0 && l < (int)height && m >= 0 && m < (int)width && l != ix && m != jy)
					{
						if(mat[l][m] != ind)
						{
							bool check = false;
							for(uint k = 0; k < boundElem.size(); ++k)
							{
								//if boundElem contains this surrounding pixel code
								if(boundElem[k] == mat[l][m])
								{
									check = true;
									boundCount[k]++;
									break;
								}
							}
							//if boundElem doesn't contain this surrounding pixel code
							if(!check)
							{
								boundElem.push_back(mat[l][m]);
								boundCount.push_back(1);
							}
						}
					}
				}
			}

			//a simple check
			if(boundCount.size() != boundElem.size())
			{
				std::cout << "Error size of boundElem and boundCount not equal" << std::endl;
				exit(EXIT_FAILURE);
			}

			//if boundElem is not empty implies a boundary pixel
			if(!boundElem.empty())
			{
				//get the boundElem with max boundCount
				uint max = 0;
				for(uint k = 1; k < boundElem.size(); ++k)
				{
					if(boundCount[k] > boundCount[max])
					{
						max = k;
					}
				}
				//dissolve pixel by changing its code in mat
				mat[ix][jy] = boundElem[max];
				Point2 temp = 
				{
					.x = (double)ix,
					.y = (double)jy
				};

				//push this pixel's coord in the assigned region
				regionPixelCoord[boundElem[max]].push_back(temp);

				//erase this pixel from current region(the one that needs to be dissolved)
				regionPixelCoord[ind].erase(regionPixelCoord[ind].begin() + i);
			}
			boundCount.clear();
			boundElem.clear();
		}
	}
}

/*
 * Returns an ImageMatrix after reverting back the whole number {0,1,..n}
 * coded bitmap to the pixel based bitmap.
 */
ImageMatrix *CodeImage::getFinalImage()
{
	ImageMatrix *m = new ImageMatrix[1];

	m->height = height;
	m->width = width;

	//initializes pixmap of ImageMatrix
	m->pixMap = new pixel*[height];
	for(uint i = 0; i < height; ++i)
	{
		m->pixMap[i] = new pixel[width];
	}

	for(uint i = 0; i < height; ++i)
	{
		for(uint j = 0; j < width; ++j)
		{
			//assign the color corresponding to code mat[i][j]
			m->pixMap[i][j] = colorCode[mat[i][j]];
		}
	}

	return m;
}