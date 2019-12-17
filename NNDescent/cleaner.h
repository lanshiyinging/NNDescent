#ifndef CLEANER_H
#define CLEANER_H

#include "nnitem.h"
#include <vector>
#include <map>

using namespace std;

class Cleaner
{
    public:
        Cleaner();
        virtual ~Cleaner();
        static void clearNbs(vector<NbHood> &nbGraph);
        static void clearKNNGraph(map<unsigned int, vector<unsigned int>* > &kNNGraph);
        static void clearKNNGraph(vector<vector<MiniNN> > &kNNGraph);

};

#endif
