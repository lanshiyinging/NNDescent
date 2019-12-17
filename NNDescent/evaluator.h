#ifndef EVALUATOR_H
#define EVALUATOR_H

/*******************************************************************
I don't want to invovle two many different blocks for one algorithm,
so I just integrate the evaluation part with this tool.

It is a pitty, this tool become complicated, which betrays my belief.

@author Wan-Lei Zhao
@date   19-Oct-2011

********************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>

using namespace std;

class Evaluator
{
    private:
        static const unsigned int dvd0;
    public:
        Evaluator(){}
        virtual ~Evaluator(){}

        static float getKGraphRecall(const char *srcFn, const char *grdFn, const unsigned int topk, const unsigned int atop);

        static int  intersect(set<unsigned int> &set1, set<unsigned int> &set2);
        static int  intersect(vector<unsigned int> &vect1, set<unsigned int> &set1, const unsigned int k);
        static int  intersect(vector<unsigned int> &vect1, set<unsigned int> &set1);
        static int  loadGrdKNNGraph(const char *srcFn, map<unsigned int, set<unsigned int>* > &kNNGraph, const int topk);
        static void loadKNNGraph(const char *srcFn, map<unsigned int, vector<unsigned int>* > &kNNGraph, const int topk, const int bd);

        static void NNRecall(const char *srcfn, const char *groundtruth,
                             const unsigned int topk0, const char *dstfn);
        static void clear_grdmap(map<string, set<string>* > &grdmap);
        static void test();
        static void testYJ();
};

#endif
