#ifndef NNDESCENT_H
#define NNDESCENT_H

#include <vector>
#include <random>

#include "nnitem.h"

/*****
This project is an implementation about NN-Descent algorithm that is presented in WWW'11 paper by Dr. Wei Dong
The implementation is fulfilled by Wan-Lei Zhao and optimized by Peng-Cheng Lin


Date: 2019-Dec-06
**/

class NNDescent
{
    private:
         unsigned int k0;
         unsigned int ndim, nCmps, ndat;
         float *data;
    private:
        vector<vector<MiniNN> > knnGraph;  ///keeps i's neighbors
        vector<NbHood> nbGraph; ///collect neighbor's and reverse neighbors for cross matching

    public:
        NNDescent();
        virtual ~NNDescent();

        int         initKnnGraph(const unsigned int k);
        int         updateLst(const unsigned int id, const unsigned int nb, float dst);
        int         appndLst(const unsigned int id,  const unsigned int nb, float dst);
        int         getRndSds(std::mt19937 &rng, const unsigned int k, unsigned int *addr, const unsigned int N);
        unsigned    nnDescent();
        int         getNBGraph(const unsigned int smplNum);
        int         buildKNNGraph(const char *srcFn, const char *dstFn, const unsigned int k);
        int         saveKNNGraph(const char *dstFn);
    public:
        static void test();
};

#endif // NNDESCENT_H
