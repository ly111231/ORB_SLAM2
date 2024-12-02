//
// Created by qingyu on 24-1-1.
//

#ifndef ORB_SLAM2_PERF_H
#define ORB_SLAM2_PERF_H

#include<iostream>
#include<algorithm>
#include<fstream>
#include<chrono>
#include<option.h>

using namespace std;

extern  int OUTPUT_DATA;

class Perf {
public:
    
    // Initialize the Perf system.
    Perf(const string& message);

//    perfReset(const string message);

//    perfReset();

    void perfStartTime();

    double perfEndTime();

    void perfSummerTime();


private:
#ifdef COMPILEDWITHC11
    std::chrono::steady_clock::time_point startTime;
#else
    std::chrono::monotonic_clock::time_point startTime;
#endif

#ifdef COMPILEDWITHC11
    std::chrono::steady_clock::time_point endTime;
#else
    std::chrono::monotonic_clock::time_point startTime;
#endif
    double maxTime,minTime,totalTime;
    int total;
    const string message;
    
};


#endif //ORB_SLAM2_PERF_H
