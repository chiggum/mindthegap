//util.h
#ifndef __UTIL_H__
#define __UTIL_H__ 1

typedef unsigned char uchar;
typedef unsigned int uint;

#include <vector>
#include <string>
#include "GraphicsGems.h" 
#include <bitset>

/*
 * Struct to store final bezier curve points corresponding to a line segment.
 * There is 1-to-1 mapping of a line segment and its corresponding curve.
 * reverse points to the reversed curve.
 * Note: pt vector stores pts in order: SC1C2EC1C2E....C1C2E
 * where S = Start point, C1 = control point 1, C2 = control point 2, E = endpoint of sub bezier curve.
 * End point becomes the start point of successor bezier curve.
 */
struct Curve
{
	Curve* reverse;
	std::vector<Point2> pt;
	Point2 start;
	Point2 end;
};

/*
 * Struct to store line segments and island line segments specs.
 * line segments = control point to control point
 * island segment = Non control point to same point 
 * path vector stores the nodes ids on the path.
 */
struct Line
{
    uint start;
    uint end;
    std::vector<uint> path;
};

/*
 * pixel (r,g,b) alpha = 1f or 255 always
 */
struct pixel
{
	uchar r;
	uchar g;
	uchar b;
};

/*
 * Struct to store region specs.
 * closedPath vector stores a vector of pointer to curve.
 * A single element of closed path contains pointers to all those curves
 * which together form a closed path.
 * curveNum stores the id = index of curves adjacent to the region.
 * col represents a pixel of the region(stores the color of the region).
 */
struct Region
{
	std::vector< std::vector<Curve*> > closedPath;
	std::vector<uint> curveNum;
	pixel col;
};

/*
 * Struct to store bitmap specs.
 * pixMap represents all the pixels of the bitmap.
 * height and width of the bitmap.
 */
struct ImageMatrix
{
	pixel **pixMap;
	uint height;
	uint width;
};
 
/*
 * Struct that represents an adjacency list node
 */
struct AdjList
{
    std::vector<uint> node;  		//ids' of nodes adjacent to this node
    uint id;						//id of node
    uint x;							//Actual x coordinate of node in image
   	uint y;							//Actual y coordinate of node in image
    std::bitset<1> isUsedUp;		//is node used in line segment formation
    std::bitset<1> isCntrlPoint;	//is node a control point
    std::bitset<1> isIslandPoint;	//is node an island point
   	std::vector<uint> adjRegion;	//Adjacent region codes
};

//returns true if pixels have same rgb values else false
bool ifEqualPixel(pixel, pixel);

//The image argument has width * height RGBA pixels or width * height * 4 bytes
void encodeOneStep(const char*, std::vector<uchar>&, uint, uint);

//Converts RGB to hexcode and return a string of hexcode
std::string RGBToHex(uint, uint, uint);

//Converts Hexcode to RGB pixel
pixel hexToRGB(std::string);

//Converts rgb(R,G,B) to RGB pixel
pixel parseRGB(std::string);

//returns true if two points (Point2) are equal else false
bool ifEqualPoint2(Point2, Point2);

#endif