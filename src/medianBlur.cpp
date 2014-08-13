//medianBlur.cpp
#include "medianBlur.h"

/*
 * Performs median blur on "m" using
 * (2 * kernelRadius + 1) X (2 * kernelRadius + 1) square kernel.
 * Final blurred image is stored back in "m".
 */
void medianBlur(int kernelRadius, ImageMatrix *m)
{
    uint h, w;
    long unsigned int L2Cache = 16 * 1024; // 16Kb: size of L2 cache

    h = m->height;
    w = m->width;
    
    uchar *src = new uchar[3*h*w];
    uchar *dst = new uchar[3*h*w];

    //flattening 2d matrix of pixels to 1d array
    for(uint i = 0; i < h; ++i)
    {
        for(uint j = 0; j < w; ++j)
        {
            src[3 * (i * w + j)+0] = m->pixMap[i][j].r;
            src[3 * (i * w + j)+1] = m->pixMap[i][j].g;
            src[3 * (i * w + j)+2] = m->pixMap[i][j].b;
        }
    }

    //call to constant time median filter
    ctmf(src, dst, w, h, 3*w, 3*w, kernelRadius, 3, L2Cache);

    //copying back results to "m"
    for(uint y = 0; y < h; y++)
    {
        for(uint x = 0; x < w; x++)
        {
            m->pixMap[y][x].r = dst[3 * (y * w + x)];
            m->pixMap[y][x].g = dst[3 * (y * w + x)+1];
            m->pixMap[y][x].b = dst[3 * (y * w + x)+2];
        }
    }

    //free memory
    delete[] src;
    delete[] dst;
}