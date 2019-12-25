cmake_minimum_required(VERSION 3.8)
project(NNDescent)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        data/sift100k.txt
        data/sift100k_gold_knn30.txt
        cleaner.cpp
        cleaner.h
        distmsr.cpp
        distmsr.h
        evaluator.cpp
        evaluator.h
        iodelegator.cpp
        iodelegator.h
        main.cpp
        Makefile
        NNDescent.cbp
        nndescent.cpp
        NNDescent.depend
        nndescent.h
        NNDescent.layout
        nnitem.h
        vstring.cpp
        vstring.h)

add_executable(NNDescent ${SOURCE_FILES})