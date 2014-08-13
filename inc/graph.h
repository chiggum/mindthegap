//graph.h
#ifndef __GRAPH_H__
#define __GRAPH_H__ 1

#include "util.h"
#include "svg.h"
#include <vector>
#include <string>

typedef Point2 *BezierCurve;

/*
 * A class of Graph using adjacency list representation
 * that adds vertices dynamically i.e. vertices (boundary pixels) are
 * added as they are discovered(see Bitmap::detectControlPoints())
 * 
 * V: Real time number of vertices (Not TOTAL until all vertices are added)
 *
 * vertex: stores all the AdjList nodes = white vertices.
 *
 * lineSeg: stores line segments from control points to control points.
 * islandLineSeg: stores line segments from NON control point to itself.
 *
 * curve: stores bezier curve(svg data) of the line segments.
 *
 * region: stores all the disjoint regions in the bitmap poppped out boundaries.
*/
class Graph
{
	private:
		uint V;	
		std::vector<AdjList> vertex;
		std::vector<Line> lineSeg;
		std::vector<Line> islandLineSeg;
		std::vector<Curve> curve;
		std::vector<Region> region;
		SVG *outputSVG;
		uint imageHeight;
		uint imageWidth;
		ImageMatrix *pop;
	public:
		/*
		 * Please refer corresponding src file graph.cpp for
		 * the functioning of the following methods.
		 */
		Graph(uint, uint, std::vector<pixel>&, ImageMatrix*);
		~Graph(); 
		void addVertex(uint, uint, uint, bool, bool, bool);
		void addEdge(uint, uint);
		void addRegion(uint, uint);
		void setControlPoint(uint, bool);
		void setIslandPoint(uint, bool);
		void formLineSegments(uint);
		void moveToNode(std::vector<uint>&, uint, uint, uint, uint, uint);
		bool checkCntrlPtAdj(uint);
		uint getAdjCntrlPtInd(uint, uint);
		bool equalSideRegions(uint, uint);
		bool isAdjToFrom(uint, uint);
		void preprocessLineSegments();
		void formCurves(double, double);
		Curve* reverseCurve(Curve*);
		void assignCurveNumToRegion();
		void processRegions();
		void assignClosedPaths(Region&);
		int getNextPathIndex(std::vector<uint>&, Point2);
		bool ifForwardDirection(int, Point2);
		std::vector<uint> DouglasPeucker(std::vector<uint>&, double);
		double shortestDistanceToSegment(uint, uint, uint);
		void writeOuput(std::string, pixel);
		bool checkIfCntrlPt(uint);
		bool checkIfUsedUp(uint);
		void removeConnection(uint, uint);
};

/*
 * Method to fit bezier curve to a series of points.
 */
void FitCurve(Point2*, int, double, double);

/*
 * Method to catch and store generated/streamed bezier curve points.
 */
void DrawBezierCurve(int, BezierCurve);

#endif