#include <iostream>

#include "nndescent.h"
#include "evaluator.h"

using namespace std;

void test()
{
    NNDescent::test();
    Evaluator::test();
}

int main()
{
    clock_t start, finish;
    double totaltime;
    start = clock();
    test();
    finish = clock();
    totaltime = (double)(finish-start)/CLOCKS_PER_SEC;
    cout << totaltime << endl;
    
    return 0;
}
