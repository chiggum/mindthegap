//svg.h
#ifndef __SVG_H__
#define __SVG_H__ 1

#include <vector>
#include <string>
#include "util.h"

/*
 * Class to draw final SVG output and output without filling(i.e only boundaries) for evaluation.
 */
class SVG
{
	private:
		uint imageHeight;
		uint imageWidth;
	public:
		/*
		 * Please refer corresponding src file svg.cpp for
		 * the functioning of the following methods.
		 */
		SVG(uint, uint);
		~SVG();
		void writeDisjointLineSegments(std::vector<Curve>&);
		void writeFinalOutput(std::vector<Region>&, std::string, pixel);
};

#endif