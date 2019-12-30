cmake_minimum_required(VERSION 3.8)
project(NNDescent)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        cmake-build-debug/cmake_install.cmake
        cmake-build-debug/CMakeCache.txt
        cmake-build-debug/Makefile
        cmake-build-debug/NNDescent
        cmake-build-debug/NNDescent.cbp
        data/sift100k.txt
        data/sift100k_gold_knn30.txt
        result/sift100k_k=40.txt
        cleaner.cpp
        cleaner.h
        CMakeLists.txt
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