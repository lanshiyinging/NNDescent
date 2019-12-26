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
    unsigned int i, j, l, s;
    for(i = 0; i < this->ndat; ++ i){
        //vector<MiniNN> nb = this->knnGraph[i];
        vector<unsigned int> newnb;
        vector<unsigned int> oldnb;
        for(j = 0; j < this->k0; ++ j){
            if(this->knnGraph[i][j].nw){
                newnb.push_back(this->knnGraph[i][j].idx);
                this->knnGraph[i][j].nw = 0;

                //reverted neighbor
                vector<MiniNN> nnb = this->knnGraph[this->knnGraph[i][j].idx];
                for(l = 0; l < this->k0; ++ l){
                    if(nnb[l].idx == i)
                        break;
                }
                if(l == this->k0)
                    this->nbGraph[this->knnGraph[i][j].idx].rnewnb.push_back(i);
            }
            else
                oldnb.push_back(this->knnGraph[i][j].idx);
            /*
            vector<MiniNN> nnb = this->knnGraph[nb[j].idx];
            for(l = 0; l < this->k0; ++ l){
                if(nnb[l].idx == i)
                    break;
            }
            if(l == this->k0)
                this->nbGraph[j].rnewnb.push_back(i);
                */
        }
        this->nbGraph[i].oldnb = oldnb;
        this->nbGraph[i].newnb = newnb;
    }
    for(i = 0; i < this->ndat; ++ i){
        int rnsize = static_cast<int>(this->nbGraph[i].rnewnb.size());
        if(this->nbGraph[i].rold == 0 && this->nbGraph[i].rnew == 0)
            continue;
        vector<unsigned int> rnewnb;
        vector<unsigned int> roldnb;
        for(j = 0; j < this->nbGraph[i].rnew; ++ j){
            this->nbGraph[i].roldnb.push_back(this->nbGraph[i].rnewnb[j]);
        }
        for(j = this->nbGraph[i].rnew; j < rnsize; ++ j){
            vector<unsigned int>::iterator iter = find(this->nbGraph[i].roldnb.begin(), this->nbGraph[i].roldnb.end(), this->nbGraph[i].rnewnb[j]);
            if(iter != this->nbGraph[i].roldnb.end())
                roldnb.push_back(this->nbGraph[i].rnewnb[j]);
            else
                rnewnb.push_back(this->nbGraph[i].rnewnb[j]);
        }
        /*
        if(rnewnb.size() + roldnb.size() > smplNum){
            std::random_device rd;
            std::mt19937 gen(rd());
            int rsize = rnewnb.size() + roldnb.size() - 1;
            std::uniform_int_distribution<> dis(0, rsize);
            vector<int> idxs;
            for(s = 0; s < smplNum; ++ s){
                int idx = dis(gen);
                while(find(idxs.begin(), idxs.end(), idx) != idxs.end())
                    idx = dis(gen);
                idxs.push_back(idx);
                if(idx < rnewnb.size())
                    this->nbGraph[i].rnewnb.push_back(rnewnb[idx]);
                else
                    this->nbGraph[i].roldnb.push_back(roldnb[idx-rnewnb.size()]);
            }
        }
         */
        this->nbGraph[i].rnewnb = rnewnb;
        this->nbGraph[i].roldnb = roldnb;
        this->nbGraph[i].rnew = static_cast<unsigned short>(rnewnb.size());
        this->nbGraph[i].rold = static_cast<unsigned short>(roldnb.size());
    }

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
    for(i = 0; i < this->ndat; ++ i){
        vector<unsigned int> newnb = this->nbGraph[i].newnb;
        vector<unsigned int> oldnb = this->nbGraph[i].oldnb;
        for(j = 0; j < this->nbGraph[i].rnew; ++ j)
            newnb.push_back(this->nbGraph[i].rnewnb[j]);
        for(j = 0; j < this->nbGraph[i].rold; ++ j)
            oldnb.push_back(this->nbGraph[i].roldnb[j]);
        //cross-comparison between old and new
        for(j = 0; j < newnb.size(); ++ j){
            for(k = 0; k < oldnb.size(); ++ k){
                dst = DistMsr::l2f(data, newnb[j], data, oldnb[k], this->ndim);
                updateLst(newnb[j], oldnb[k], dst);
                updateLst(oldnb[k], newnb[j], dst);
                ccmps += 1;
            }
        }
        //cross-comparison in-between new
        for(a = 0; a < newnb.size(); ++ a){
            for(b = a+1; b < newnb.size(); ++ b){
                dst = DistMsr::l2f(data, newnb[a], data, newnb[b], this->ndim);
                updateLst(newnb[a], newnb[b], dst);
                updateLst(newnb[b], newnb[a], dst);
                ccmps += 1;
            }
        }

    }

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
        this->nCmps += ccmps;
        Cleaner::clearNbs(this->nbGraph);
        i++;
    }while(i < 6 && ccmps > 512);
    cout<< i << "\n";
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
    const char *srcFn = "/home/sylan/lanshiying/multimedia/project2/NNDescent/NNDescent/data/sift100k.txt";
    const char *dstFn = "/home/sylan/lanshiying/multimedia/project2/NNDescent/NNDescent/result/sift100k_k=40.txt";

    NNDescent *mynn = new NNDescent();
    int k = 40;
    mynn->buildKNNGraph(srcFn, dstFn, k);

}
