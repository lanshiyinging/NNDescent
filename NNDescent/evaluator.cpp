#include "evaluator.h"

#include "iodelegator.h"
#include "cleaner.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <cmath>
#include <map>

/**
*
*
Evaluate quality of k-NN graph
*
*
**/

using namespace std;

const unsigned int Evaluator::dvd0 = 10000;

float Evaluator::getKGraphRecall(const char *srcFn, const char *grdFn, const unsigned int topk, const unsigned int atop)
{
    float r = 0;
    unsigned int key = 0, counts = 0, nhits = 0, hit = 0, sz = 0;
    map<unsigned int, vector<unsigned int>* > kNNGraph;
    map<unsigned int, set<unsigned int>* > grdNNGraph;
    map<unsigned int, vector<unsigned int>* >::iterator mit;
    vector<unsigned int>* crntVect;
    set<unsigned int>* crntSet;
    set<unsigned int>::iterator sit;
    vector<unsigned int>::iterator vit;

    sz = Evaluator::loadGrdKNNGraph(grdFn, grdNNGraph, atop);
    Evaluator::loadKNNGraph(srcFn, kNNGraph, atop, sz);

    cout<<"Size is "<<" ................... "<<sz<<" done\n";
    for(mit = kNNGraph.begin(); mit != kNNGraph.end(); mit++)
    {
        key = mit->first;
        crntVect = mit->second;
        crntSet  = grdNNGraph[key];
        ///cout<<key<<"\t"<<crntSet->size()<<endl;
        if(crntSet != NULL)
        {
             hit = Evaluator::intersect(*crntVect, *crntSet, atop);
             ///if(ri == 1 && key <= 100)
             ///cout<<key<<"\t"<<hit<<endl;
             nhits += hit;
             counts++;
        }
    }
    r  = nhits/(counts*atop+0.0);
    cout<<"Recall is: "<<r<<"\tnHits: "<<nhits<<endl;

    //Cleaner::clearKNNGraph(grdNNGraph);
    Cleaner::clearKNNGraph(kNNGraph);

    return r;
}

int Evaluator::loadGrdKNNGraph(const char *srcFn, map<unsigned int, set<unsigned int>* > &kNNGraph, const int topk)
{
    unsigned int id = 0, cid = 0;
    unsigned int sz, it, dim = 0;
    int i = 0, numb = 0;
    set<unsigned int>* crntItms = NULL;

    kNNGraph.clear();
    ifstream *inStrm = new ifstream(srcFn, ios::in);
    if(!inStrm->is_open())
    {
        cout<<"File '"<<srcFn<<"' cannot open for read!\n";
        exit(0);
    }

    (*inStrm)>>sz;
    (*inStrm)>>dim;
    for(it = 0; it < sz; it++)
    {
        (*inStrm)>>id;
        (*inStrm)>>numb;
        crntItms = new set<unsigned int>;
        kNNGraph.insert(pair<unsigned int, set<unsigned int>*>(id, crntItms));
        for(i = 0; i < numb; i++)
        {
            (*inStrm)>>cid;
            if(i < topk)
            crntItms->insert(cid);
        }
    }

    inStrm->close();

    return sz;
}


void Evaluator::loadKNNGraph(const char *srcFn, map<unsigned int, vector<unsigned int>* > &kNNGraph, const int topk, const int bd)
{
    int sz = 0, it = 0, i = 0, numb = 0;
    unsigned int id = 0, cid = 0;
    vector<unsigned int>* crntItms = NULL;

    kNNGraph.clear();
    ifstream *inStrm = new ifstream(srcFn, ios::in);
    if(!inStrm->is_open())
    {
        cout<<"File '"<<srcFn<<"' cannot open for read!\n";
        exit(0);
    }

    (*inStrm)>>sz;
    (*inStrm)>>numb;

    sz = sz < bd?sz:bd;
    for(it = 0; it < sz; it++)
    {
        (*inStrm)>>id;
        (*inStrm)>>numb;

        crntItms = new vector<unsigned int>;
        kNNGraph.insert(pair<unsigned int, vector<unsigned int>*>(id, crntItms));
        for(i = 0; i < numb; i++)
        {
            (*inStrm)>>cid;
            if(i < topk)
            crntItms->push_back(cid);
        }
    }

    inStrm->close();

    return ;
}

int Evaluator::intersect(set<unsigned int> &set1, set<unsigned int> &set2)
{
    set<unsigned int>::iterator sit;
    int cmmn = 0;
    for(sit = set1.begin(); sit != set1.end(); sit++)
    {
        if(set2.find(*sit) != set2.end())
        {
            cmmn++;
        }
    }
    return cmmn;
}

int Evaluator::intersect(vector<unsigned int> &vect1, set<unsigned int> &set1, const unsigned int k)
{
     int ovrlap = 0;
     unsigned int i = 0;
     for(i = 0; i < k; i++)
     {
         if(set1.find(vect1[i]) != set1.end())
         {
              ovrlap++;
         }
     }
     return ovrlap;
}

int Evaluator::intersect(vector<unsigned int> &vect1, set<unsigned int> &set1)
{
     int ovrlap = 0;
     vector<unsigned int>::iterator vit;
     for(vit = vect1.begin(); vit != vect1.end(); vit++)
     {
         if(set1.find(*vit) != set1.end())
         {
              ovrlap++;
         }
     }
     return ovrlap;
}

void Evaluator::clear_grdmap(map<string, set<string>* > &grdmap)
{
    map<string, set<string> *>::iterator mit;
    set<string>::iterator sit;
    set<string> *crnt_lst;

    for(mit = grdmap.begin(); mit != grdmap.end(); mit++)
    {
        string key = mit->first;
        crnt_lst   = mit->second;
        for(sit = crnt_lst->begin(); sit != crnt_lst->end(); sit++)
        {
            string mkey = *sit;
            mkey.clear();
        }
        key.clear();
        crnt_lst->clear();
        delete crnt_lst;
    }
    grdmap.clear();
    return ;
}


void Evaluator::test()
{
    const char *knngrd1  = "/home/sylan/lanshiying/multimedia/project2/NNDescent/NNDescent/data/sift100k_gold_knn30.txt";
    const char *knnsrc1  = "/home/sylan/lanshiying/multimedia/project2/NNDescent/NNDescent/result/sift100k_k=40.txt";
    const char *knngrd2  = "/home/wlzhao/datasets/gist1m_gold_knn30.txt";
    const char *knnsrc2  = "/home/wlzhao/datasets/gist1m_k=40.txt";
    const char *knngrd3  = "/home/wlzhao/datasets/rand1m_gold_knn30.txt";
    const char *knnsrc3  = "/home/wlzhao/datasets/rand1m_k=40.txt";

    const unsigned int topk2 = 10;

    Evaluator::getKGraphRecall(knnsrc1, knngrd1, 1, 1);
    Evaluator::getKGraphRecall(knnsrc1, knngrd1, 1, topk2);

    return ;
}

