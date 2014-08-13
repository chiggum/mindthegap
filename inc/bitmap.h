//bitmap.h
#ifndef __BITMAP_H__
#define __BITMAP_H__ 1

#include "util.h"
#include "codeImage.h"
#include "graph.h"
#include "lodepng.h"
#include <vector>
#include <string>
#include "info.h"

/*
 * Class to store all bitmaps assosciated with an image.
 * orig: original bitmap
 * 
 * pre: preprocessed bitmap which is same as the original image with a 
 * row at top and bottom and a col at left and right of black pixels(0,0,0)
 * 
 * pop: popped out boundaries bitmap is an image with boundaries of different 
 * regions popped out and represented as white pixels (255,255,255)
 * 
 * codedImage stores the pointer to an object of codeImage which holds the corresponding
 * coded map of pop image.
 * 
 * graph represents a Graph formed by all white vertices of pop.
 * 
 * pixToNodeMap stores the map of a pixel to the id of Node in graph. If a pixel is
 * not white then it maps to -1 and dangerous pixels are mapped to -2.
 * Please refer graph::formLineSegments for the definition of dangerous points.
 *
 * boundaryPixel: pixel used in popped out boundaries.
 *
 * borderPixel: pixel used in appended rows and cols.
 *
 * inputParam: holds the input parameter(option) values provided by the user.
 * 
 * Note: pop image has height and width twice that of pre image (See Algorithm 
 * for popping out boundary for more details).
 */
class Bitmap
{
	private:
		ImageMatrix* orig;
		ImageMatrix* pre;
		ImageMatrix* pop;
		CodeImage *codedImage;
		Graph *graph;
		int **pixToNodeMap;
		uint numControlPoints;
		pixel boundaryPixel;
		pixel borderPixel;
		Info *inputParam;
	public:
		/*
		 * Please refer corresponding src file bitmap.cpp for
		 * the functioning of the following methods.
		 */
		Bitmap(class Info*);
		~Bitmap();
		void del(ImageMatrix*);
		void decodeOneStep(const char*, uint&, uint&, std::vector<uchar>&);
		void intialize(ImageMatrix**, uint, uint);
		void preprocess();
		void popoutBoundaries();
		void detectControlPoints();
		int checkUniqueRegionPixel(pixel, std::vector<pixel>&);
		void formAdjacencyList();
		void processImage();
		void writeOuputSVG();
		void removeNoise();
		void removeDangerousConnections();
};

#endif