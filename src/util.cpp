//util.cpp
#include "util.h"
#include "lodepng.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cctype>

//Returns true if pixel x and pixel y have same rgb values else false
bool ifEqualPixel(pixel x, pixel y)
{
	if(x.r == y.r && x.g == y.g && x.b == y.b)
		return true;
	else
		return false;
}

//The image argument has width * height RGBA pixels or width * height * 4 bytes
void encodeOneStep(const char* filename, std::vector<uchar>& image, unsigned width, unsigned height)
{
  //Encode the image
  unsigned error = lodepng::encode(filename, image, width, height);

  //if there's an error, display it
  if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
}

//Input: r,g,b values
//Return corresponding hexcode
std::string RGBToHex(uint rNum, uint gNum, uint bNum)
{
    std::string result;
    //start with #
    result.append("#");

    char buff[16];
    sprintf(buff, "%02X%02X%02X",rNum,gNum,bNum);
    result.append(buff);
    return result;
}

//Input: hash/hex coded (base 16) rgb value in form #RRGGBB
//Return corresponding rgb pixel
pixel hexToRGB(std::string hexCode)
{
    std::string hex = "";
    pixel ret;
    std::string legit = "0123456789abcdefABCDEF";

    //will always be true becuase it is prechecked
    if(hexCode.substr(0, 1) == "#")
    {
      hex = hexCode.substr(1);
    }

    //check if hex contatins legible hex characters
    for(uint i = 0; i < hex.size(); ++i)
    {
      //if hex[i] is not found among legible characters
      if(std::string::npos == legit.find_first_of(hex[i]))
      {
        std::cout << "Please check your hex color code characters" << std::endl;
        exit(EXIT_FAILURE);
      }
    }

    std::string r = "", g = "", b = "";
    r = hex.substr(0, 2);
    g = hex.substr(2, 2);
    b = hex.substr(4, 2);

    ret.r = (uchar)strtol(r.c_str(), NULL, 16);
    ret.g = (uchar)strtol(g.c_str(), NULL, 16);
    ret.b = (uchar)strtol(b.c_str(), NULL, 16);

    return ret;
}

//Input: base 10 rgb value in form rgb(R,G,B)
//Return corresponding rgb pixel
pixel parseRGB(std::string rgb)
{
    std::string trimRGB = "";
    char *pEnd;
    pixel ret;

    //will always be true because it is prechecked
    if(rgb.substr(0, 3) == "rgb")
    {
      trimRGB = rgb.substr(3);
    }

    //replace non numeric characters(excluding negative sign '-') with blank space
    for(uint i = 0; i < trimRGB.size(); ++i)
    {
      if(!isdigit(trimRGB[i]) && trimRGB[i] != '-')
      {
        trimRGB[i] = ' ';
      }
    }

    uint r = strtol(trimRGB.c_str(), &pEnd, 10);
    uint g = strtol(pEnd, &pEnd, 10);
    uint b = strtol(pEnd, NULL, 10);

    //if r/g/b value doesn't lie b/w [0, 255]
    if(!(r >= 0 && r <= 255 && g >= 0 && g <= 255 && b >= 0 && b <= 255))
    {
      std::cout << "Please check your R G B values: must lie in [0, 255]" << std::endl;
      exit(1);
    }
    ret.r = (uchar)r;
    ret.g = (uchar)g;
    ret.b = (uchar)b;

    return ret;
}

//Returns true if a and b have same (x, y) values else false
bool ifEqualPoint2(Point2 a, Point2 b)
{
  if(a.x == b.x && a.y == b.y)
    return true;
  else
    return false;
}