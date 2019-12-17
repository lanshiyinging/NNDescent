#include "distmsr.h"

DistMsr::DistMsr()
{
    //ctor
}

DistMsr::~DistMsr()
{
    //dtor
}

float DistMsr::l2f(const float v1[], const unsigned int s1,
                  const float v2[], const unsigned int s2, const unsigned int d0)
{
    float w1 = 0.0f, w2 = 0.0f;
    int loc1 = d0 * s1;
    int loc2 = d0 * s2;

    for(unsigned int i = 0; i < d0; i++)
    {
        w1  = v1[loc1+i] - v2[loc2+i];
        w2 += w1*w1;
    }

    ///w2 = sqrt(w2);

    return w2;
}
