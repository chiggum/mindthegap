//medianBlur.h
#ifndef __MEDIANBLUR_H__
#define __MEDIANBLUR_H__ 1

#include "util.h"
#include "ctmf.h"

/*
 * Performs median blur on "m" using
 * (2 * kernelRadius + 1) X (2 * kernelRadius + 1) square kernel.
 * Final blurred image is stored back in "m".
 */
void medianBlur(int kernelRadius, ImageMatrix *m);

#endif