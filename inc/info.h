//info.h
#ifndef __INFO_H__
#define __INFO_H__ 1

#include <string>
#include "util.h"

/*
 * Class to store the options value(info) and parse arguments
 * to main provided by the user.
 *
 * toleranceCurve: tolerance used by the bezier curve fitting code.
 * It is the max tolerated distance between a coordinate on the original line segment
 * (i.e. series of points) and its corresponding mapped coordinate on the fitted bezier curve.
 *
 * toleranceLine: max tolerated shortest distance from the coordinate on
 * the original line segment (i.e. series of points) to the fitted straight line.
 *
 * turdSize: Those regions which have size <= turdSize are removed from the bitmap
 * while preprocessing noisy image. size of a region = no. of pixels in that region.
 *
 * switchNoisy: true implies noisy image processing will be done.
 *
 * medianBlurKernelSize(n):square kernel of size nXn is used for median blur
 * of noisy image.
 *
 * KMeansMaxIters: maximum no. of iterations used in K Means clustering algorithm 
 * while performing posterization of median blurred image.
 *
 * numClusters: Number of clusters to be made in K Means clustering algorithm.
 */
class Info
{
	public:
		/*
		 * All members are public
		 */
		std::string inputFileName;
		std::string outFileName;
		std::string color;
		pixel bgColor;
		double toleranceCurve;
		double toleranceLine;
		bool bgColorProvided;
		bool switchNoisy;
		int turdSize;
		int medianBlurKernelSize;
		int KMeansMaxIters;
		int numClusters;

		/*
		 * Please refer corresponding src file info.cpp for
		 * the functioning of the following methods.
		 */
		Info();
		~Info();
		void parseInputArg(int, char**);
};

#endif