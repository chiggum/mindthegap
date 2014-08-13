//graph.cpp
#include "graph.h"
#include "debug.h"
#include <fstream>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <bitset>
#include <string>
#include "util.h"


std::vector<Point2> ptStore; //stores bezier points as they are genereted/streamed

//LOOP: length of the loop from a control point to itself to make it to be considered as line segment
#define LOOP 10
//LIMIT: length of the segment to be considered for either DouglasPeuckar simplification or jump simplification
#define LIMIT 10
//JUMP: Number of pixels to be left while walking over the line segment
#define JUMP 0
//GAP: Net absolute change in the x and y direction of two pixels to be considered for a closed loop
#define GAP 150
//EPSILON: For douglas pecker algorithm
#define EPSILON 1

/*
 * Method which is called whenever a bezier curve is generated.
 * A bezier curve will have 4 pts. S, C1, C2, E.
 * To catch all Bezier curves assosciated with a line segment, 
 * ptStore global variable is used.
 * n = number of points generated ALWAYS 3
 */
void DrawBezierCurve(int n, BezierCurve curve)
{
  // only the first three points are stored (S C1 C2),
  // the end point is the start of the next curve
  for(int i = 0; i < n; ++i) {
    Point2 temp = {
      .x = curve[i].x,
      .y = curve[i].y
    };
    ptStore.push_back(temp);
  }
}

//Constructor
Graph::Graph(uint h, uint w, std::vector<pixel> &v, ImageMatrix *x)
{
	V = 0;
	imageHeight = h;
	imageWidth = w;
	outputSVG = new SVG(imageHeight, imageWidth);

	//initializing region with the col value set.
	for(uint i = 0; i < v.size(); ++i)
	{
		Region *tempRegion = new Region[1];
		tempRegion->col = v[i];
		region.push_back(*tempRegion);
	}

	pop = x;
}

//Destructor
Graph::~Graph()
{
	//delete vertex
	for(uint i = 0; i < vertex.size(); ++i)
	{
		vertex[i].node.clear();
		vertex[i].adjRegion.clear();
	}
	vertex.clear();

	//delete lineSeg
	for(uint i = 0; i < lineSeg.size(); ++i)
	{
		lineSeg[i].path.clear();
	}
	lineSeg.clear();

	//delete islandLineSeg
	for(uint i = 0; i < islandLineSeg.size(); ++i)
	{
		islandLineSeg[i].path.clear();
	}
	islandLineSeg.clear();

	//delete outputSVG
	delete outputSVG;
}

//Adds a vertex to the graph with specified parameters set
void Graph::addVertex(uint id, uint x, uint y, bool isUsedUp = 0, bool isCntrlPoint = 0, bool isIslandPoint = 0)
{
	AdjList *list = new AdjList[1];
	list->id = id;
	list->x = x;
	list->y = y;
	list->isUsedUp.set(0,isUsedUp);
	list->isCntrlPoint.set(0,isCntrlPoint);
	list->isIslandPoint.set(0,isIslandPoint);
	vertex.push_back(*list);
	V = V + 1;	
}

// Adds an edge to graph
void Graph::addEdge(uint src, uint dest)
{
    // Add an edge from src to dest.  A new node is added to the adjacency list of src.
    vertex[src].node.push_back(dest);
}

//Adds the code of adjacent region(ind) to the adjacent region vector of source vertex(src)
void Graph::addRegion(uint src, uint ind)
{
	vertex[src].adjRegion.push_back(ind);
}

//Sets isControlPoint flag of vertex[a]
void Graph::setControlPoint(uint a, bool b)
{
	vertex[a].isCntrlPoint.set(0,b);
}

//Sets isIslandPoint flag of vertex[a]
void Graph::setIslandPoint(uint a, bool b)
{
	vertex[a].isIslandPoint.set(0,b);
}


/*
 * 	formLineSegments and MoveToNode
 *
 * 	It forms line segments from control point to control point and island line segments too.
 * 	Pass = 1 forms linesegments
 * 	Pass = 2 forms islandLineSegments
 *	Note: Before forming line segments, "dangerous connections" surrounding dangerous points must be removed.
 *
 *	Dangerous points/pixels: those non boundary pixels which have boundary pixels at UP, DOWN, LEFT, RIGHT.
 *	Dangerous connections:
 *		Let x be a dangerous point. 
 *		If, in popped out boundary bitmap, x->UPLEFT and x->DOWNRIGHT have same
 *		rgb value as x then connection between UP and LEFT & DOWN and RIGHT are dangerous.
 *		Else if, in popped out boundary bitmap, x->UPRIGHT and x->DOWNLEFT have same
 *		rgb value as x then connection between UP and RIGHT & DOWN and LEFT are dangerous.
 *		Following are the pictorial reresentation of above two cases:
 *			x 	B 	?			? 	B 	x
 *			B 	x 	B 			B 	x 	B
 *			? 	B 	x 			x 	B 	?
 *	
 *	A high level view of the algorithm.
 *
 * 	pass 1: forms line segments from control point to control point.
 * 		start from a control point(X) and initialize an empty "path" variable
 * 		move away via adjacencly list represnetation of graph and a visited node/vertex bookkeeping
 * 		once you touch another control point OR (same control point with length of path > some FIXED CONSTANT)
 * 			push all the visited nodes/vertices from the starting cp to this cp to the "path"
 * 			push this "path" in a vector of lineSeg
 * 		repeat from step 2 with the same starting control point(X) but different adjacent node/vertex which is not yet visited
 * 		repeat the above steps until all control points are exhausted
 *
 * 		Now, you have a vector of lineSeg which contains all line segments from cp to cp.
 *
 * 	pass 2: forms island line segments from an unused node/vertex to itself
 * 		start from an unsed node/vertex and initialize an empty "path" variable
 * 		move away via adjacencly list represnetation of graph and a visited node/vertex bookkeeping
 * 		once you touch the same starting node/vertex
 * 			push all the visited nodes/vertices from starting pt to itself to the "path"
 * 			push this "path" in a vector of islandLineSeg
 * 		repeat from step 2 with the same starting point(X) but different adjacent node/vertex which is not yet visited(rarely happens)
 * 		repeat the above steps until all vertices/nodes are used up
 *
 * 		Now, you have a vector of islandLineSeg which contains all line segments from a non control point to itself.
 */
void Graph::formLineSegments(uint pass)
{
	for(uint i = 0; i < V; ++i)
	{
		if((vertex[i].isCntrlPoint.test(0) && pass == 1) || (!vertex[i].isUsedUp.test(0) && pass == 2))
		{
			if(pass == 1)
				vertex[i].isUsedUp.set(0, true);	//set usedUp true
			else
				vertex[i].isIslandPoint.set(0, true);	//set islandPoint true

			for(uint j = 0; j < vertex[i].node.size(); ++j)
			{
				//if adjacentt vertex not used up
				if(!vertex[vertex[i].node[j]].isUsedUp.test(0))
				{
					//initialize an empty line "temp"
					Line *temp = new Line[1];

					//setting this adjacent node used up i.e. visited
					vertex[vertex[i].node[j]].isUsedUp.set(0, true);

					//move or step over this adjacent node
					moveToNode(temp->path, vertex[i].node[j], i, i, pass, 1);

					//if the path is not empty
					if(!temp->path.empty())
					{
						temp->path.push_back(vertex[i].node[j]);
						temp->path.push_back(i);

						//reorienting path (reversed paths are returned)
						std::reverse(temp->path.begin(), temp->path.end());

						temp->start = temp->path.front();
						temp->end = temp->path.back();

						if(pass == 1)
							lineSeg.push_back(*temp);
						else 
						{
							islandLineSeg.push_back(*temp);
							vertex[i].isUsedUp.set(0, true);
						}
					}
					else if(pass == 2)
					{
						//if no island line segment is formed then make this adjacent node not used up
						vertex[vertex[i].node[j]].isUsedUp.set(0, false);
					}
				}
			}
		}
	}
	
/******************************************************************/
/***********************TEST 5:DEBUGGING MODE**********************/
/******************************************************************/

#ifdef _TEST_5_

	std::string path = ROOT_DIR;
	std::ofstream ofsTest5((path + "/check/test_5_cpp.txt").c_str(), std::ofstream::out);
	ofsTest5 << lineSeg.size() + islandLineSeg.size() << std::endl;

	ofsTest5 << "Line segments: " << std::endl;

	for(uint i = 0; i < lineSeg.size(); ++i)
	{
		ofsTest5 << lineSeg[i].path.size() << std::endl;
		for(uint j = 0; j < lineSeg[i].path.size(); ++j)
		{
			ofsTest5 << lineSeg[i].path[j] << " ";
		}
		ofsTest5 << std::endl;
	}

	ofsTest5 << "IsLand Line segments: " << std::endl;

	for(uint i = 0; i < islandLineSeg.size(); ++i)
	{
		ofsTest5 << islandLineSeg[i].path.size() << std::endl;
		for(uint j = 0; j < islandLineSeg[i].path.size(); ++j)
		{
			ofsTest5 << islandLineSeg[i].path[j] << " " ;
		}
		ofsTest5 << std::endl;
	}

	ofsTest5.close();

#endif

/******************************************************************/
/***********************TEST 5:DEBUGGING MODE END******************/
/******************************************************************/

/******************************************************************/
/***********************EVAL 3:DEBUGGING MODE**********************/
/******************************************************************/

#ifdef _EVAL_3_

	#ifndef _TEST_5_
		std::string path = ROOT_DIR;
	#endif
	//generate some image
	std::vector<unsigned char> image;
	image.resize(imageWidth * imageHeight * 4);
	for(uint y = 0; y < imageHeight; y++)
	{
		for(uint x = 0; x < imageWidth; x++)
		{
			image[4 * imageWidth * y + 4 * x + 0] = pop->pixMap[y][x].r;
			image[4 * imageWidth * y + 4 * x + 1] = pop->pixMap[y][x].g;
			image[4 * imageWidth * y + 4 * x + 2] = pop->pixMap[y][x].b;
			image[4 * imageWidth * y + 4 * x + 3] = 255;
		}
	}

	for(uint i = 0; i < lineSeg.size(); ++i)
	{
		for(uint j = 0; j < lineSeg[i].path.size(); ++j)
		{
			image[4 * imageWidth * vertex[lineSeg[i].path[j]].x + 4 * vertex[lineSeg[i].path[j]].y + 0] = 0;
			image[4 * imageWidth * vertex[lineSeg[i].path[j]].x + 4 * vertex[lineSeg[i].path[j]].y + 1] = 0;
			image[4 * imageWidth * vertex[lineSeg[i].path[j]].x + 4 * vertex[lineSeg[i].path[j]].y + 2] = 0;
			image[4 * imageWidth * vertex[lineSeg[i].path[j]].x + 4 * vertex[lineSeg[i].path[j]].y + 3] = 255;			
		}
	}

	for(uint i = 0; i < islandLineSeg.size(); ++i)
	{
		for(uint j = 0; j < islandLineSeg[i].path.size(); ++j)
		{
			image[4 * imageWidth * vertex[islandLineSeg[i].path[j]].x + 4 * vertex[islandLineSeg[i].path[j]].y + 0] = 255;
			image[4 * imageWidth * vertex[islandLineSeg[i].path[j]].x + 4 * vertex[islandLineSeg[i].path[j]].y + 1] = 0;
			image[4 * imageWidth * vertex[islandLineSeg[i].path[j]].x + 4 * vertex[islandLineSeg[i].path[j]].y + 2] = 0;
			image[4 * imageWidth * vertex[islandLineSeg[i].path[j]].x + 4 * vertex[islandLineSeg[i].path[j]].y + 3] = 255;			
		}
	}

	for(uint i = 0; i < V; ++i)
	{
		if(!vertex[i].isUsedUp.test(0))
		{
			image[4 * imageWidth * vertex[i].x + 4 * vertex[i].y + 0] = 255;
			image[4 * imageWidth * vertex[i].x + 4 * vertex[i].y + 1] = 255;
			image[4 * imageWidth * vertex[i].x + 4 * vertex[i].y + 2] = 0;
			image[4 * imageWidth * vertex[i].x + 4 * vertex[i].y + 3] = 255;
		}
	}

	encodeOneStep((path + "/check/eval_3_lineSeg_beforeprocess.png").c_str(), image, imageWidth, imageHeight);
	image.clear();

#endif

/******************************************************************/
/***********************EVAL 3:DEBUGGING MODE END******************/
/******************************************************************/

}

/*
 * Moves to a node.
 * See the comments before formLineSegments
 */
void Graph::moveToNode(std::vector<uint> &v, uint to, uint from, uint startedFrom, uint pass, uint len)
{
	uint ind;

	if(pass == 2 && to == startedFrom)
	{
		v.push_back(to);
	}
	else if(pass == 1 && checkCntrlPtAdj(to) && ((!vertex[ind = getAdjCntrlPtInd(to, startedFrom)].isUsedUp.test(0)) || (len > LOOP && ind == startedFrom)))
	{
		v.push_back(ind);
	}
	else
	{
		for(uint i = 0; i < vertex[to].node.size(); ++i)
		{
			//if adjacent regions are equal && not used up / not visited && not adjacent to "from" vertex && not equal to from
			if(equalSideRegions(to, vertex[to].node[i]) && !vertex[vertex[to].node[i]].isUsedUp.test(0) && !isAdjToFrom(from, vertex[to].node[i]) && vertex[to].node[i] != from)
			{
				std::vector<uint> temp;

				//set used up / visited true
				vertex[vertex[to].node[i]].isUsedUp.set(0, true);

				//increment length of path
				len++;

				//move to adjacent node
				moveToNode(temp, vertex[to].node[i], to, startedFrom, pass, len);

				//if temporary path "temp" is not empty
				if(!temp.empty())
				{
					v.assign(temp.begin(), temp.end());

					if(vertex[to].node[i] != startedFrom)
						v.push_back(vertex[to].node[i]);
					//else no need to add because starting point already pushed in
					
					//no need to move to adjacent node is temp has control point in the beginning i.e. at End after reorientation
					if(vertex[v[0]].isCntrlPoint.test(0))
						break;
				}
				else if(pass == 2)
				{
					//make vertex unused if an island line segment is not formed
					vertex[vertex[to].node[i]].isUsedUp.set(0, false);
				}
			}
		}
	}
}

/*
 * Returns true if there is an adjacent control point to vertex[a] else false
 */
bool Graph::checkCntrlPtAdj(uint a)
{
	for(uint i = 0; i < vertex[a].node.size(); ++i)
	{
		if(vertex[vertex[a].node[i]].isCntrlPoint.test(0))
			return true;
	}
	return false;
}

/*
 * Returns the index of the adjacent control point
 * if more than two adjacent control points then returns index other than startedfrom
 * vertex[startedfrom] is always a control point here
 */
uint Graph::getAdjCntrlPtInd(uint to, uint startedFrom)
{
	std::vector<uint> v;
	for(uint i = 0; i < vertex[to].node.size(); ++i)
	{
		if(vertex[vertex[to].node[i]].isCntrlPoint.test(0))
			v.push_back(vertex[to].node[i]);
	}
	if(v.size() == 1)
		return v[0];
	else
	{
		for(uint i = 0; i < v.size(); ++i)
		{
			if(v[i] != startedFrom)
				return v[i];
		}
	}
	//Will never come at this point because it is already checked that there is an adjacent control point
	std::cout << "Error while retrieving adjacent control point: Graph::getAdjCntrlPtInd" << std::endl;
	exit(1);
}

/*
 * Returns true if the regions adjacent to vertex[a] and vertex[b]
 * are equal else false.
 */
bool Graph::equalSideRegions(uint a, uint b)
{
	if(vertex[a].adjRegion.size() != vertex[b].adjRegion.size())
		return false;

	//sort the region codes
	std::sort(vertex[a].adjRegion.begin(), vertex[a].adjRegion.end());
	std::sort(vertex[b].adjRegion.begin(), vertex[b].adjRegion.end());

	//compare the sorted region codes
	for(uint i = 0; i < vertex[a].adjRegion.size(); ++i)
	{
		if(vertex[a].adjRegion[i] != vertex[b].adjRegion[i])
			return false;
	}

	return true;
}

/*
 * returns true if vertex[dest] is adjacent to vertex[from] else false
 */
bool Graph::isAdjToFrom(uint from, uint dest)
{
	for(uint i = 0; i < vertex[from].node.size(); ++i)
	{
		if(vertex[from].node[i] == dest)
			return true;
	}
	return false;
}

/*
 * Preprocesses line segments with either taking points at equal intervals on a line segment
 * OR using Douglas Peucker algorithm
 */
void Graph::preprocessLineSegments()
{
	
	//Preprocess : Take points at an appropriate interval depending on the length of line segment OR DouglasPeucker
#if 1
	if(JUMP != 0)
	{
		std::vector<uint> tempLineSeg;
		for(uint i = 0; i < lineSeg.size(); ++i)
		{
			tempLineSeg.clear();
			if(lineSeg[i].path.size() > LIMIT)// if line length is above 6 then only go in and preprocess 
			{
				tempLineSeg.push_back(lineSeg[i].path[0]);
				tempLineSeg.push_back(lineSeg[i].path[1]);
				for(uint j = 2; j < lineSeg[i].path.size() - 2; j = j + 1 + JUMP)
				{
					tempLineSeg.push_back(lineSeg[i].path[j]);
				}
				tempLineSeg.push_back(lineSeg[i].path[lineSeg[i].path.size() - 2]);
				tempLineSeg.push_back(lineSeg[i].path[lineSeg[i].path.size() - 1]);

				lineSeg[i].path.clear();
				lineSeg[i].path.assign(tempLineSeg.begin(), tempLineSeg.end());
			}
		}

		for(uint i = 0; i < islandLineSeg.size(); ++i)
		{
			tempLineSeg.clear();
			if(islandLineSeg[i].path.size() > LIMIT)// if line length is above 5 then only go in and preprocess 
			{
				tempLineSeg.push_back(islandLineSeg[i].path[0]);
				tempLineSeg.push_back(islandLineSeg[i].path[1]);
				for(uint j = 2; j < islandLineSeg[i].path.size() - 2; j = j + 1 + JUMP)
				{
					tempLineSeg.push_back(islandLineSeg[i].path[j]);
				}
				tempLineSeg.push_back(islandLineSeg[i].path[islandLineSeg[i].path.size() - 2]);
				tempLineSeg.push_back(islandLineSeg[i].path[islandLineSeg[i].path.size() - 1]);

				islandLineSeg[i].path.clear();
				islandLineSeg[i].path.assign(tempLineSeg.begin(), tempLineSeg.end());
			}
		}

		tempLineSeg.clear();
	}
#endif
#if 0

	std::vector<uint> tempLineSeg;
	for(uint i = 0; i < lineSeg.size(); ++i)
	{
		tempLineSeg.clear();
		if(lineSeg[i].path.size() > LIMIT)// if line length is above 5 then only go in and preprocess 
		{
			tempLineSeg = DouglasPeucker(lineSeg[i].path,EPSILON);

			lineSeg[i].path.clear();
			lineSeg[i].path.push_back(lineSeg[i].start);
			lineSeg[i].path.insert(lineSeg[i].path.begin() + 1, tempLineSeg.begin(), tempLineSeg.end());
			lineSeg[i].path.push_back(lineSeg[i].end);
		}
	}

	for(uint i = 0; i < islandLineSeg.size(); ++i)
	{
		tempLineSeg.clear();
		if(islandLineSeg[i].path.size() > LIMIT)// if line length is above 5 then only go in and preprocess 
		{
			tempLineSeg.push_back(islandLineSeg[i].path[0]);
			tempLineSeg.push_back(islandLineSeg[i].path[1]);
			for(uint j = 2; j < islandLineSeg[i].path.size() - 2; j = j + 1 + JUMP)
			{
				tempLineSeg.push_back(islandLineSeg[i].path[j]);
			}
			tempLineSeg.push_back(islandLineSeg[i].path[islandLineSeg[i].path.size() - 2]);
			tempLineSeg.push_back(islandLineSeg[i].path[islandLineSeg[i].path.size() - 1]);

			islandLineSeg[i].path.clear();
			islandLineSeg[i].path.assign(tempLineSeg.begin(), tempLineSeg.end());
		}
	}

	tempLineSeg.clear();
#endif


/******************************************************************/
/***********************EVAL 3:DEBUGGING MODE**********************/
/******************************************************************/

#ifdef _EVAL_3_

	//generate some image
	std::string path = ROOT_DIR;
	std::vector<unsigned char> image;
	image.resize(imageWidth * imageHeight * 4);
	for(uint y = 0; y < imageHeight; y++)
	{
		for(uint x = 0; x < imageWidth; x++)
		{
			image[4 * imageWidth * y + 4 * x + 0] = pop->pixMap[y][x].r;
			image[4 * imageWidth * y + 4 * x + 1] = pop->pixMap[y][x].g;
			image[4 * imageWidth * y + 4 * x + 2] = pop->pixMap[y][x].b;
			image[4 * imageWidth * y + 4 * x + 3] = 255;
		}
	}

	for(uint i = 0; i < lineSeg.size(); ++i)
	{
		for(uint j = 0; j < lineSeg[i].path.size(); ++j)
		{
			image[4 * imageWidth * vertex[lineSeg[i].path[j]].x + 4 * vertex[lineSeg[i].path[j]].y + 0] = 0;
			image[4 * imageWidth * vertex[lineSeg[i].path[j]].x + 4 * vertex[lineSeg[i].path[j]].y + 1] = 0;
			image[4 * imageWidth * vertex[lineSeg[i].path[j]].x + 4 * vertex[lineSeg[i].path[j]].y + 2] = 0;
			image[4 * imageWidth * vertex[lineSeg[i].path[j]].x + 4 * vertex[lineSeg[i].path[j]].y + 3] = 255;			
		}
	}

	for(uint i = 0; i < islandLineSeg.size(); ++i)
	{
		for(uint j = 0; j < islandLineSeg[i].path.size(); ++j)
		{
			image[4 * imageWidth * vertex[islandLineSeg[i].path[j]].x + 4 * vertex[islandLineSeg[i].path[j]].y + 0] = 255;
			image[4 * imageWidth * vertex[islandLineSeg[i].path[j]].x + 4 * vertex[islandLineSeg[i].path[j]].y + 1] = 0;
			image[4 * imageWidth * vertex[islandLineSeg[i].path[j]].x + 4 * vertex[islandLineSeg[i].path[j]].y + 2] = 0;
			image[4 * imageWidth * vertex[islandLineSeg[i].path[j]].x + 4 * vertex[islandLineSeg[i].path[j]].y + 3] = 255;			
		}
	}

	encodeOneStep((path + "/check/eval_3_lineSeg_afterprocess.png").c_str(), image, imageWidth, imageHeight);
	image.clear();

#endif

/******************************************************************/
/***********************EVAL 3:DEBUGGING MODE END******************/
/******************************************************************/

}

std::vector<uint> Graph::DouglasPeucker(std::vector<uint> &v, double epsilon)
{
	double dmax = 0;
	uint index = 0;
	for(uint i = 2; i < v.size() - 1; ++i)
	{
		double d = shortestDistanceToSegment(v[i], v.front(), v.back());
		if(d > dmax)
		{
			index = i;
			dmax = d;
		}
	}

	std::vector<uint> v3;

	if(dmax > epsilon)
	{
		std::vector<uint> v1(v.begin(), v.begin() + index);
		std::vector<uint> v2(v.begin() + index, v.end());
		v3 = DouglasPeucker(v1, epsilon);
		std::vector<uint> v4 = DouglasPeucker(v2, epsilon);
		v3.insert(v3.end(), v4.begin(), v4.end());
	}
	else
	{
		v3.push_back(v.front());
		v3.push_back(v.back());
	}
	return v3;
}

double Graph::shortestDistanceToSegment(uint i, uint j, uint k)
{
	double x1 = vertex[i].x;
	double y1 = vertex[i].y;
	double x2 = vertex[j].x;
	double y2 = vertex[j].y;
	double x3 = vertex[k].x;
	double y3 = vertex[k].y;
	if(x2 == x3)
		return abs(x1 - x2);

	double m = (y3 - y2) / (x3 - x2);
	double c = y2 - m * x2;

	return abs((y1 - m * x1 - c) / (sqrt(1 + m * m)));
}

/*
 * Form curves by fitting line segments to bezier curve.
 * The coordinates are rescaled before sending them to fitting routine.
 */
void Graph::formCurves(double toleranceCurve, double toleranceLine)
{
	for(uint i = 0; i < lineSeg.size(); ++i)
	{
		Point2 *d = new Point2[lineSeg[i].path.size()];

		for(uint j = 0; j < lineSeg[i].path.size(); ++j)
		{
			d[j].x = (vertex[lineSeg[i].path[j]].y-2)*0.5;	//rescaling
			d[j].y = (vertex[lineSeg[i].path[j]].x-2)*0.5;	//rescaling
		}

		ptStore.clear();
		Curve *tempCurve = new Curve[1];
		tempCurve->reverse = new Curve[1];

		//call fitting routine on the series of points in d
		FitCurve(d, lineSeg[i].path.size(), toleranceCurve, toleranceLine);

		ptStore.push_back(d[lineSeg[i].path.size()-1]);
		tempCurve->pt = ptStore;	// this is a copy by value
		tempCurve->start = ptStore.front();
		tempCurve->end = ptStore.back();
		tempCurve->reverse = reverseCurve(tempCurve);
		curve.push_back(*tempCurve);
	}

	for(uint i = 0; i < islandLineSeg.size(); ++i)
	{
		Point2 *d = new Point2[islandLineSeg[i].path.size()];

		for(uint j = 0; j < islandLineSeg[i].path.size(); ++j)
		{
			d[j].x = (vertex[islandLineSeg[i].path[j]].y-2)*0.5;	//rescaling
			d[j].y = (vertex[islandLineSeg[i].path[j]].x-2)*0.5;	//rescaling
		}


		ptStore.clear();
		Curve *tempCurve = new Curve[1];
		tempCurve->reverse = new Curve[1];

		//call fitting routine on the series of points in d
		FitCurve(d, islandLineSeg[i].path.size(), toleranceCurve, toleranceLine);

		ptStore.push_back(d[islandLineSeg[i].path.size()-1]);
		tempCurve->pt = ptStore; // this is a copy by value
		tempCurve->start = ptStore.front();
		tempCurve->end = ptStore.back();
		tempCurve->reverse = reverseCurve(tempCurve);
		curve.push_back(*tempCurve);
	}


/******************************************************************/
/***********************TEST 6:DEBUGGING MODE**********************/
/******************************************************************/

#ifdef _TEST_6_

	std::string path = ROOT_DIR;
	std::ofstream ofsTest6((path + "/check/test_6_cpp.txt").c_str(), std::ofstream::out);
	ofsTest6 << curve.size() << std::endl;

	for(uint i = 0; i < curve.size(); ++i)
	{
		ofsTest6 << curve[i].pt.size() << std::endl;
		ofsTest6 << curve[i].start.x << " " << curve[i].start.y << std::endl;
		ofsTest6 << curve[i].end.x << " " << curve[i].end.y << std::endl;
		for(uint j = 0; j < curve[i].pt.size(); ++j)
		{
			ofsTest6 << curve[i].pt[j].x << " " << curve[i].pt[j].y << std::endl;
		}
	}

	ofsTest6.close();

#endif

/******************************************************************/
/***********************TEST 6:DEBUGGING MODE END******************/
/******************************************************************/


/******************************************************************/
/***********************EVAL 4:DEBUGGING MODE**********************/
/******************************************************************/

#ifdef _EVAL_4_

	outputSVG->writeDisjointLineSegments(curve);

#endif

/******************************************************************/
/***********************EVAL 4:DEBUGGING MODE END******************/
/******************************************************************/
}

/*
 * Reverses a curve
 */
Curve* Graph::reverseCurve(Curve *x)
{
	Curve *p =  new Curve[1];
	p->reverse = x;

	p->start = x->end;
	p->end = x->start;

	//Never true
	if(x->pt.size() < 3)
		std::cout << "Error: Curve has less than 3 control points: Graph::reverseCurve" << std::endl;

	for(int i = x->pt.size()-1; i >= 0; --i)
	{
		p->pt.push_back(x->pt[i]);
	}

	return p;
}

/*
 * Assigns the indices of those curves to a region which are adjacent to that region
 */
void Graph::assignCurveNumToRegion()
{
	for(uint i = 0; i < lineSeg.size(); ++i)
	{
		for(uint j = 0; j < lineSeg[i].path.size(); ++j)
		{
			//if lineSeg[i].path[j] is not a control point
			if(!vertex[lineSeg[i].path[j]].isCntrlPoint.test(0))
			{
				//if there are two adj Regions to this point
				if(vertex[lineSeg[i].path[j]].adjRegion.size() == 2)
				{
					region[vertex[lineSeg[i].path[j]].adjRegion[0]].curveNum.push_back(i);
					region[vertex[lineSeg[i].path[j]].adjRegion[1]].curveNum.push_back(i);
					break;
				}
			}
		}
	}
	for(uint i = 0; i < islandLineSeg.size(); ++i)
	{
		for(uint j = 0; j < islandLineSeg[i].path.size(); ++j)
		{
			//if there are two adj Regions to this point
			if(vertex[islandLineSeg[i].path[j]].adjRegion.size() == 2)
			{
				region[vertex[islandLineSeg[i].path[j]].adjRegion[0]].curveNum.push_back(i + lineSeg.size());
				region[vertex[islandLineSeg[i].path[j]].adjRegion[1]].curveNum.push_back(i + lineSeg.size());
				break;
			}
		}
	}


/******************************************************************/
/***********************TEST 7:DEBUGGING MODE**********************/
/******************************************************************/

#ifdef _TEST_7_

	std::string path = ROOT_DIR;
	std::ofstream ofsTest7((path + "/check/test_7_cpp.txt").c_str(), std::ofstream::out);
	ofsTest7 << "No. of regions" << region.size() << std::endl;

	for(uint i = 0; i < region.size(); ++i)
	{
		ofsTest7 << "No. of curves assigned to region " << i << " : " << region[i].curveNum.size() << std::endl;
		for(uint j = 0; j < region[i].curveNum.size(); ++j)
		{
			ofsTest7 << region[i].curveNum[j] << " ";
		}
		ofsTest7 << std::endl;
	}

	ofsTest7.close();
	
#endif

/******************************************************************/
/***********************TEST 7:DEBUGGING MODE END******************/
/******************************************************************/

}

/*
 * processes all regions
 */
void Graph::processRegions()
{
	for(uint i = 0; i < region.size(); ++i)
	{
		assignClosedPaths(region[i]);
	}
}

/*
 * write output svg to outFileName by transferring control to writeFinalOutput method of SVG class.
 */
void Graph::writeOuput(std::string outFileName, pixel bg)
{
	outputSVG->writeFinalOutput(region, outFileName, bg);
}

/*
 * Assigns paths surrounding a region in a proper order which ultimately forms a closed path taken together.
 *
 * Basic idea: 
 * We have a pool of line segments that surrounds a region.
 * We need to arrange them in proper order so that, taken in order, they form a closed path.
 * This is achieved by first picking a random line segment.
 * Then the one which has starting point same/nearest to end point of last line segment added, is attached further.
 * Now, above step is repeated until all line segments are used up.
 */
void Graph::assignClosedPaths(Region &rgn)
{
	bool firstTime = true;
	int focus = 0;
	std::vector<Curve*> *tempPath = NULL;
	bool initiallyEmpty = rgn.curveNum.empty();

	//until all surrounding paths are consumed
	while(!rgn.curveNum.empty())
	{
		if(firstTime)
		{
			//take any initial path (say the first/last one from curveNum)
			tempPath = new std::vector<Curve*>();
			tempPath->clear();
			focus = rgn.curveNum.back();
			rgn.curveNum.pop_back();
			tempPath->push_back(&curve[focus]);
			firstTime = false;
		}
		else
		{
			//if start and point of last path added are same i.e. island
			if(tempPath->front()->start.x == tempPath->back()->end.x && tempPath->front()->start.y == tempPath->back()->end.y)
			{
				rgn.closedPath.push_back(*tempPath);
				firstTime = true;
			}
			else
			{
				//get next path/line which can be joined with the end point of last path/line added
				focus = getNextPathIndex(rgn.curveNum, tempPath->back()->end);
				if(focus == -1)
					std::cout << "Error: focus = - 1 i.e. rgn.curveNum is empty still in the loop : Graph::assignClosedPaths" << std::endl;
				else
				{
					//Erasing the corresponding curve number
					uint temp = rgn.curveNum[focus];
					rgn.curveNum.erase(rgn.curveNum.begin() + focus);
					focus = temp;

					//if forward direction path has to be attached
					if(ifForwardDirection(focus, tempPath->back()->end))
					{
						tempPath->push_back(&curve[focus]);
					}
					else
					{
						tempPath->push_back(curve[focus].reverse);
					}
				}
			}
		}
	}
	//if curveNum is not empty initially then push the ordered paths in the closedpath of rgn	
	if(!initiallyEmpty)
		rgn.closedPath.push_back(*tempPath);
}

/*
 * Get next path index whose start or end point is nearest to the point b
 */
int Graph::getNextPathIndex(std::vector<uint> &v, Point2 b)
{
	double val = 1000;  //No distance will exceed this disatance
	int ind = -1;

	for(uint i = 0; i < v.size(); ++i)
	{
		double myX = abs(curve[v[i]].start.x - b.x);
		double myY = abs(curve[v[i]].start.y - b.y);
		if(myX + myY < val)
		{
			val = myX + myY;
			ind = i;
		}

		myX = abs(curve[v[i]].end.x - b.x);
		myY = abs(curve[v[i]].end.y - b.y);
		if(myX + myY < val)
		{
			val = myX + myY;
			ind = i;
		}
	}
	
	return ind;
}

/*
 * Checks if curve's start point is near to "b" or end point
 * if start point then return true else false
 */
bool Graph::ifForwardDirection(int focus, Point2 b)
{
	double val1;
	double val2;
	val1 = abs(curve[focus].start.x - b.x) + abs(curve[focus].start.y - b.y);
	val2 = abs(curve[focus].end.x - b.x) + abs(curve[focus].end.y - b.y);

	if(val1 < val2)
		return true;
	else
		return false;
}

//Returns true if vertex[ind] is a control point else false
bool Graph::checkIfCntrlPt(uint ind)
{
	return vertex[ind].isCntrlPoint.test(0);
}

//Returns true if vertex[ind] is used up else false
bool Graph::checkIfUsedUp(uint ind)
{
	return vertex[ind].isUsedUp.test(0);
}

//Remove connection/edge between vertex[a] and vertex[b]
void Graph::removeConnection(uint a, uint b)
{
	vertex[a].node.erase(std::remove(vertex[a].node.begin(), vertex[a].node.end(), b), vertex[a].node.end());
	vertex[b].node.erase(std::remove(vertex[b].node.begin(), vertex[b].node.end(), a), vertex[b].node.end());
}