#ifndef DISTMSR_H
#define DISTMSR_H


class DistMsr
{
    public:
        DistMsr();
        virtual ~DistMsr();

    static float l2f(const float v1[], const unsigned int s1,
              const float v2[], const unsigned int s2, const unsigned int d0);
};

#endif // DISTMSR_H
