//codeImage.h
#ifndef __CODEIMAGE_H__
#define __CODEIMAGE_H__ 1

#include "util.h"
#include "GraphicsGems.h"
#include <vector>

/*
 * Classs to create and store the coded image matrix
 * i.e. assigning different unsigned integers to different
 * disjoint regions of image.
 *
 * mat: stores image coded with ints from {0, 1, .. n}; n = number of regions.
 *
 * colorCode: stores the rgb value value corresponding to the region/code.
 *
 * regionPixelCoord[i]: vector containing coordinates of pixels in ith region. 
 * regionPixelCoord: vector contatining above vectors corresponding to each region.
 */
class CodeImage
{
	private:
		int** mat;	
		uint height;
		uint width;
		uint numDisjointRegions;
		std::vector<pixel> colorCode;
		std::vector< std::vector<Point2> > regionPixelCoord;	
	public:
		/*
		 * Please refer corresponding src file codeImage.cpp for
		 * the functioning of the following methods.
		 */
		CodeImage(ImageMatrix*);
		~CodeImage();
		void codeImage(ImageMatrix*);
		int **getMatrix();
		std::vector<pixel>* getColCode();
		void processRegions(uint);
		void dissolveRegion(int);
		ImageMatrix *getFinalImage();
};

#endif