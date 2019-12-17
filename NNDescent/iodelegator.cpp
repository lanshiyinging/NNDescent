#include <iostream>
#include <cstring>
#include <cassert>
#include <iomanip>
#include <cstdio>
#include <cmath>

#include "iodelegator.h"
#include "vstring.h"
//#include "vmath.h"

using namespace std;

long IODelegator::getFileSize(const char *srcfn)
{
    ifstream *inStrm = new ifstream(srcfn, ios::in|ios::binary);
    inStrm->seekg(0, ios::end);
    long sz = (long)inStrm->tellg();
    inStrm->close();
    delete inStrm;

    return sz;
}

float *IODelegator::load_refSet(const char *srcFn, unsigned int &d, unsigned int &r)
{
    float *refMat = NULL;

    if(VString::endWith(srcFn, ".fvecs"))
    {
        refMat = IODelegator::load_fvecs(srcFn, d, r);
    }
    else if(VString::endWith(srcFn, ".txt"))
    {
        refMat = IODelegator::loadMatrix(srcFn, r, d);
        ///cout<<r<<"\t"<<d<<endl;
    }else
    {
         cout<<"Unrecorgnizable file format!\n";
    }
    return refMat;
}

float  *IODelegator::loadMatrix(const char *fn, unsigned int &row, unsigned int &col)
{
    unsigned int irow = 0, idim = 0, loc = 0;
    float vals[2] = {0};
    assert(fn);

    ifstream *inStrm = new ifstream();
    inStrm->open(fn, ios::in);
    if(inStrm->fail())
    {
        cout<<"Fail to read "<<fn<<endl;
        delete inStrm;
        exit(0);
    }

    (*inStrm)>>vals[0];
    (*inStrm)>>vals[1];
    //cout<<vals[0]<<"\t"<<vals[1]<<endl;
    row = (int)round(vals[0]);
    col = (int)round(vals[1]);
    float *mat = new float[row*col];


    for(irow = 0; irow < row; irow++)
    {
        loc = irow*col;
        for(idim = 0; idim < col; idim++)
        {
            (*inStrm) >>mat[loc + idim];
        }
    }
    inStrm->close();
    delete inStrm;
    return mat;
}

void IODelegator::saveMatrix(const char *dstFn, unsigned int row, unsigned int col, const float *mat)
{
    unsigned int irow = 0, idim = 0, loc = 0;
    assert(dstFn);

    ofstream *outStrm = new ofstream();
    outStrm->open(dstFn, ios::out);
    if(outStrm->fail())
    {
        cout<<"Fail to read '"<< dstFn << "'" <<endl;
        delete outStrm;
        exit(0);
    }

    (*outStrm)<<row<<" ";
    (*outStrm)<<col<<endl;

    for(irow = 0; irow < row; irow++)
    {
        loc = irow*col;
        (*outStrm) << mat[loc];
        for(idim = 1; idim < col; idim++)
        {
            (*outStrm) << " " << mat[loc + idim];
        }
        (*outStrm)<<endl;
    }
    outStrm->close();
    delete outStrm;
    return ;
}


float *IODelegator::load_fvecs(const char *fvecFn, unsigned int &d, unsigned int &r)
{
    float *mat = NULL, *vect = NULL, *ppmat = NULL;
    unsigned bfsz = 0, bg = 0;
    unsigned int di = 0, line = 0;

    ifstream *inStrm = new ifstream(fvecFn, ios::in|ios::binary);
    if(!inStrm->is_open())
    {
        cout<<"File '"<<fvecFn<<"' cannot open for read!\n";
        return NULL;
    }
    r = d = 0;
    bg = inStrm->tellg();
    inStrm->read((char*)&di,  sizeof(unsigned int));
    d = di;
    bfsz = d*sizeof(float);
    inStrm->seekg(0, ios::end);
    r = ((long)inStrm->tellg() - bg)/(bfsz + sizeof(unsigned int));
    inStrm->close();

    mat   = new float[r*d];
    vect  = new float[d];
    ppmat = mat;

    inStrm = new ifstream(fvecFn, ios::in|ios::binary);
    while(!inStrm->eof() && line < r)
    {
        inStrm->read((char*)&di, sizeof(int));
        assert(di == d);

        bfsz = d*sizeof(float);
        inStrm->read((char*)vect, bfsz);

        line++;
        memcpy(ppmat, vect, bfsz);
        ppmat = ppmat + d;
    }

    delete [] vect;
    inStrm->close();

    return mat;
}


void IODelegator::test()
{

}


