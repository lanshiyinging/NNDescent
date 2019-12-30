#include "cleaner.h"

Cleaner::Cleaner()
{
    //ctor
}

Cleaner::~Cleaner()
{
    //dtor
}

void Cleaner::clearNbs(vector<NbHood> &nbGraph)
{
    unsigned int i = 0;

    for(i = 0; i < nbGraph.size(); i++)
    {
        nbGraph[i].newnb.clear();
        nbGraph[i].oldnb.clear();
        nbGraph[i].rnewnb.clear();
        nbGraph[i].roldnb.clear();
        nbGraph[i].rold = 0;
        nbGraph[i].rnew = 0;
    }
}

void Cleaner::clearKNNGraph(map<unsigned int, vector<unsigned int>* > &kNNGraph)
{
    map<unsigned int, vector<unsigned int>* >::iterator mit;
    vector<unsigned int>* crntItms = NULL;

    for(mit = kNNGraph.begin(); mit != kNNGraph.end(); mit++)
    {
        crntItms = mit->second;
        if(crntItms != NULL)
        {
          crntItms->clear();
          delete crntItms;
        }
    }
    kNNGraph.clear();
}

void Cleaner::clearKNNGraph(vector<vector<MiniNN> > &kNNGraph)
{
    vector<vector<MiniNN> >::iterator vit;

    for(vit = kNNGraph.begin(); vit != kNNGraph.end(); vit++)
    {
        vector<MiniNN> &crntItms = *vit;
        crntItms.clear();
    }
    kNNGraph.clear();
}
