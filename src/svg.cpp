//svg.cpp
#include "svg.h"
#include <fstream>
#include <string>
#include "util.h"

//Constructor
SVG::SVG(uint h, uint w)
{
	imageHeight = (h-4)*0.5;
	imageWidth = (w-4)*0.5;
}

//Destructor
SVG::~SVG()
{

}

/*
 * Writes an svg with curves only.
 * Both forward and reverse curve are written and the resultant
 * svg is expected with only red curves which overlap black curves
 * NOTE: This method will be called only if _EVAL_4_ is defined
 */
void SVG::writeDisjointLineSegments(std::vector<Curve> &v)
{
	std::string path = ROOT_DIR;
	std::ofstream ofsTest6((path + "/check/eval_4_borders.svg").c_str(), std::ofstream::out);

	ofsTest6 << "<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" height=\"" << imageHeight << "\" width=\"" << imageWidth << "\">" << std::endl << "<g>" << std::endl;
	//white backgoround
	ofsTest6 << "<rect width=\"" << imageWidth << "\" height=\"" << imageHeight << "\" fill=\"#ffffff\" stroke-width=\"0\" /> " << std::endl;

	//iterate over all curves
	for(uint i = 0; i < v.size(); ++i)
	{
		ofsTest6 << "<path fill=\"none\" stroke=\"black\" stroke-width=\"0.5\" d =\" ";

		//print forward curve path
		for(uint j = 0;j < v[i].pt.size(); ++j)
		{
			if(j == 0)
			{
				ofsTest6 << "M" << v[i].start.x << " " << v[i].start.y << " C";
			}
			else
				ofsTest6 << v[i].pt[j].x << "," << v[i].pt[j].y << " ";
		}
		ofsTest6 << " \" />" << std::endl;
	}

	//iterate over all curves
	for(uint i = 0; i < v.size(); ++i)
	{
		ofsTest6 << "<path fill=\"none\" stroke=\"red\" stroke-width=\"0.5\" d =\" ";

		//print reverse curve path
		for(uint j = 0;j < v[i].reverse->pt.size(); ++j)
		{
			if(j == 0)
			{
				ofsTest6 << "M" << v[i].reverse->start.x << " " << v[i].reverse->start.y << " C";
			}
			else
				ofsTest6 << v[i].reverse->pt[j].x << "," << v[i].reverse->pt[j].y << " ";
		}
		ofsTest6 << " \" />" << std::endl;
	}

	ofsTest6 << "</g>" << std::endl << "</svg>" << std::endl;

	ofsTest6.close();
}

/*
 * Writes Final SVG corresponding to original image.
 * Arguments: Vector of regions in the image, output File name and background color.
 */
void SVG::writeFinalOutput(std::vector<Region> &rgn, std::string outFileName, pixel bgColor)
{
	std::ofstream ofsFinal(outFileName.c_str(), std::ofstream::out);

	//header
	ofsFinal << "<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" height=\"" << imageHeight << "\" width=\"" << imageWidth << "\">" << std::endl << "<g>" << std::endl;
	
	//background
	ofsFinal << "<rect width=\"" << imageWidth << "\" height=\"" << imageHeight << "\" fill=\"" << RGBToHex((uint)bgColor.r, (uint)bgColor.g, (uint)bgColor.b) << "\"  stroke-width=\"0\" /> " << std::endl;

	//iterate over all regions
	for(uint i = 0; i < rgn.size(); ++i)
	{
		//No printing of empty paths
		if(rgn[i].closedPath.empty())
			continue;
		
		if(rgn[i].closedPath.size() == 1)
			ofsFinal << "<path  fill=\"" << RGBToHex((uint)rgn[i].col.r, (uint)rgn[i].col.g, (uint)rgn[i].col.b) << "\" stroke-width=\"0\" stroke=\"" << RGBToHex((uint)rgn[i].col.r, (uint)rgn[i].col.g, (uint)rgn[i].col.b) << "\" d = \"";
		else
			ofsFinal << "<path  fill-rule=\"evenodd\" fill=\"" << RGBToHex((uint)rgn[i].col.r, (uint)rgn[i].col.g, (uint)rgn[i].col.b) << "\" stroke-width=\"0\" stroke=\"" << RGBToHex((uint)rgn[i].col.r, (uint)rgn[i].col.g, (uint)rgn[i].col.b) << "\" d = \"";	

		//iterate over the closedpaths surrounding the region
		for(uint j = 0; j < rgn[i].closedPath.size(); ++j)
		{
			//iterate over the curves in a closed path
			for(uint k = 0; k < rgn[i].closedPath[j].size(); ++k)
			{
				int count = 0;				//count == 0 signifies start of a bezier curve segment
				bool cMayBePlaced = false;

				//iterate over the points in the curve
				for(uint m = 0; m < rgn[i].closedPath[j][k]->pt.size(); ++m)
				{	
					count = count % 3;					
					if(m == 0)
					{
						//print starting point of first curve in a closed path
						if(k == 0)
						{
							//if SC1C2E forms a straight line and its the start of the bezier curve segment i.e. count==0
							if(count == 0 && m+4 <= rgn[i].closedPath[j][k]->pt.size() && ifEqualPoint2(rgn[i].closedPath[j][k]->pt[m], rgn[i].closedPath[j][k]->pt[m+1]) && ifEqualPoint2(rgn[i].closedPath[j][k]->pt[m+2], rgn[i].closedPath[j][k]->pt[m+3]))
							{
								ofsFinal << " M" << rgn[i].closedPath[j][k]->start.x << " " << rgn[i].closedPath[j][k]->start.y << " L";
								m = m + 2;				//jump to the C2 skipping C1
								cMayBePlaced = true;	//C may be placed after E
							}
							else
							{
								ofsFinal << " M" << rgn[i].closedPath[j][k]->start.x << " " << rgn[i].closedPath[j][k]->start.y << " C";
								count = count + 1;
								cMayBePlaced = false;	//C is already placed
							}
						}
					}
					else
					{
						//if SC1C2E forms a straight line and its a start of the bezier curve segment i.e. count==0
						if(count == 0 && m+4 <= rgn[i].closedPath[j][k]->pt.size() && ifEqualPoint2(rgn[i].closedPath[j][k]->pt[m], rgn[i].closedPath[j][k]->pt[m+1]) && ifEqualPoint2(rgn[i].closedPath[j][k]->pt[m+2], rgn[i].closedPath[j][k]->pt[m+3]))
						{
							ofsFinal << rgn[i].closedPath[j][k]->pt[m].x << "," << rgn[i].closedPath[j][k]->pt[m].y << " L";
							m = m + 2;				//jump to the C2 skipping C1
							cMayBePlaced = true;	//C may be placed after E
						}
						else
						{
							//we don't want to place C after last bezier curve segment
							if(cMayBePlaced && !(k == rgn[i].closedPath[j].size() - 1 && m == rgn[i].closedPath[j][k]->pt.size() - 1))
								ofsFinal << rgn[i].closedPath[j][k]->pt[m].x << "," << rgn[i].closedPath[j][k]->pt[m].y << " C";
							else
								ofsFinal << rgn[i].closedPath[j][k]->pt[m].x << "," << rgn[i].closedPath[j][k]->pt[m].y << " ";
							count = count + 1;
							cMayBePlaced = false;	//C already placed or last bezier curve segment
						}
					}
				}
			}
		}
		//for closed path
		ofsFinal << " Z\" />" << std::endl;
		ofsFinal << std::endl;
	}
	//footer
	ofsFinal << "</g>" << std::endl << "</svg>" << std::endl;

	ofsFinal.close();
}