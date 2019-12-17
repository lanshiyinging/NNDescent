#ifndef IODELEGATOR_H
#define IODELEGATOR_H

#include <unordered_map>
#include <fstream>
#include <vector>
#include <string>
#include <list>
#include <map>
#include <set>

//#include "invtitem.h"

using namespace std;

/*****************************************************
In charge of input operations

@author:     Wan-Lei Zhao
@updated:    Feb-2015
@institute:  XMU.CN

*****************************************************/

class IODelegator
{
private:
    static const unsigned int LONGCHAR  = 2000;
    static const unsigned int FNLEN     = 1024;
    static const unsigned int precision = 3;

public:
    IODelegator();
    static long   getFileSize(const char *srcfn);
    static float *loadMatrix(const char *fn, unsigned int &num, unsigned int &dim);
    static void   saveMatrix(const char *dstFn, unsigned int row, unsigned int col, const float *mat);

    static float *load_refSet(const char *srcFn,   unsigned int &d, unsigned int &r);
    static float *load_fvecs(const char *fvecFn,   unsigned int &d, unsigned int &r);


    static void test();

    ~IODelegator();
};

#endif
