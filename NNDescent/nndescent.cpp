#include <iostream>
#include <fstream>
#include <cassert>
#include <ctime>

#include "iodelegator.h"
#include "nndescent.h"
#include "distmsr.h"
#include "cleaner.h"


NNDescent::NNDescent()
{
    this->data = NULL;
    this->ndat = this->ndim = 0;
}

int NNDescent::getRndSds(std::mt19937 &rng, const unsigned int k, unsigned int *addr, const unsigned int N)
{
    unsigned int i = 0;
    if (N == k)
    {
        for (i = 0; i < k; ++i)
        {
            addr[i] = i;
        }
        return 0;
    }
    for (i = 0; i < k; ++i)
    {
        addr[i] = rng() % (N - k);
    }
    sort(addr, addr + k);
    for (i = 1; i < k; ++i)
    {
        if (addr[i] <= addr[i-1])
        {
            addr[i] = addr[i-1] + 1;
        }
    }
    unsigned off = rng() % N;
    for (i = 0; i < k; ++i)
    {
        addr[i] = (addr[i] + off) % N;
    }

    return 0;
}

int NNDescent::initKnnGraph(const unsigned int k)
{
    unsigned i = 0, j = 0, L = 0;
    unsigned *seeds  = new unsigned[512];
    std::mt19937 rng(time(NULL));
    ///cout<<"Initialize nn graph .............. ";
    for( i = 0; i < ndat; i++)
    {
        knnGraph[i].resize(k);
        this->getRndSds(rng, k+1, seeds, ndat);
        L = 0;
        for (j = 0; j < k+1;  j++)
        {
            MiniNN &nb = knnGraph[i][L];

            if(seeds[j] == i) continue;

            nb.idx  = seeds[j];
            assert(nb.idx != i);
            nb.val  = DistMsr::l2f(data, i, data, nb.idx, this->ndim);
            nb.nw   = 1;
            this->nCmps++;
            L++;
            if(L >= k)
                break;
        }
        sort(knnGraph[i].begin(), knnGraph[i].end(), MiniNN::LLcomparer);
        ///cout<<i<<"\t"<<knnGraph[i].size()<<endl;
    }
    delete [] seeds;
    seeds = NULL;
    cout<<"done \n";

    return 0;
}

int NNDescent::getNBGraph(const unsigned int smplNum)
{
    ///Collect the old neighbors and new neighbors from each sample's neighborhood (including reverse neighbors)
    ///put them into nbGraph

    return 0;
}

int NNDescent::updateLst(const unsigned int id, const unsigned int nb, float dst)
{
    int i =0, topk =0, j = 0, l;
    i = topk = k0;
    vector<MiniNN> & addr = this->knnGraph[id];

    if(addr[topk-1].val < dst)
        return 0;

    while(i > 0)
    {
        j = i-1;
        if(addr[j].val <= dst) break;
        i = j;
    }
    l = i;
    while(l > 0)
    {
        j = l - 1;
        if(addr[j].val < dst) break;
        if(addr[j].idx == nb) return 0;
        l = j;
    }

    j = topk - 1;

    while(j > i)
    {
        addr[j].idx = addr[j-1].idx;
        addr[j].val = addr[j-1].val;
        addr[j].nw  = addr[j-1].nw;
        --j;
    }

    addr[j].idx = nb;
    addr[j].val = dst;
    addr[j].nw  = 1;

    return j;
}

unsigned NNDescent::nnDescent()
{
    unsigned int i = 0, j = 0, k = 0, a = 0, b = 0;
    unsigned ccmps = 0;
    float dst = 0;

    ///Perform cross-comparison between old and new neighbors, and in-between new neighbors
    ///Call updateLst(.) to insert newly produced <a, b, dst> tuples into knnGraph

    return ccmps;
}

int NNDescent::buildKNNGraph(const char *srcFn, const char *dstFn, const unsigned int k)
{
    unsigned int i = 0, ccmps = 0;
    float rate  = 0.0f;
    this->nCmps = 0;
    this->data  = IODelegator::load_refSet(srcFn, this->ndim, this->ndat);
    this->k0    = k;

    if(this->ndim == 0 || this->ndat == 0)
    exit(0);

    this->knnGraph.resize(this->ndat);
    this->initKnnGraph(k);
    this->nbGraph.resize(this->ndat);
    cout<<"start NN-Descent ................. ";
    do
    {
        getNBGraph(100);
        ccmps = nnDescent();
        Cleaner::clearNbs(this->nbGraph);
        i++;
    }while(i < 6 && ccmps > 512);
    cout<<"\n";
    rate = (nCmps*2.0)/(ndat*(ndat-1));
    cout<<"The scanning rate is: "<<rate<<endl;

    this->saveKNNGraph(dstFn);

    return 0;
}

int NNDescent::saveKNNGraph(const char *dstFn)
{
    ofstream *outStrm = new ofstream(dstFn, ios::out);
    unsigned int i = 0, j =  0;
    unsigned int nb = 20 > k0?k0:20;
    (*outStrm)<<this->knnGraph.size()<<" "<<nb<<endl;
    for(i = 0; i < this->knnGraph.size(); i++)
    {
        (*outStrm)<<i<<" "<<nb;
        for(j = 0; j < nb; j++)
        {
            (*outStrm)<<" "<<knnGraph[i][j].idx;
        }
        (*outStrm)<<endl;
    }
    outStrm->close();
    return 0;
}

NNDescent::~NNDescent()
{
    delete [] this->data;
    this->data = NULL;
    Cleaner::clearKNNGraph(this->knnGraph);
    Cleaner::clearNbs(this->nbGraph);
}

void NNDescent::test()
{
    const char *srcFn = "/home/wlzhao/datasets/bignn/mirproj/sift100k.txt";
    const char *dstFn = "/home/wlzhao/datasets/bignn/mirproj/hi.txt";

    NNDescent *mynn = new NNDescent();
    int k = 40;
    mynn->buildKNNGraph(srcFn, dstFn, k);

}
