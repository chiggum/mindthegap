//ctmf.h
#ifndef __CTMF_H__
#define __CTMF_H__ 1

//constant time median filter
//please refer ctmf.cpp for more detail
void ctmf(
        const unsigned char* const, unsigned char* const,
        const int, const int,
        const int, const int,
        const int, const int, const long unsigned int
        );

#endif