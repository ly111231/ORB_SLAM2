//
// Created by qingyu on 24-1-1.
//
#include "Perf.h"
    int OUTPUT_DATA = 0;
    Perf::Perf(const string& message):maxTime(-1),minTime(-1),totalTime(0),total(0),message(message){};

    void Perf::perfStartTime(){
        #ifdef COMPILEDWITHC11
            startTime = std::chrono::steady_clock::now();
        #else
            startTime = std::chrono::monotonic_clock::now();
        #endif
    }

    double Perf::perfEndTime(){
        #ifdef COMPILEDWITHC11
            endTime = std::chrono::steady_clock::now();
        #else
            endTime = std::chrono::monotonic_clock::now();
        #endif

        double curTime = std::chrono::duration_cast<std::chrono::duration<double> >(endTime - startTime).count();

        if(maxTime <= 0 || maxTime < curTime){
            maxTime = curTime;
        }

        if(minTime <= 0 || minTime > curTime){
            minTime = curTime;
        }

        total++;
        totalTime += curTime;

        return curTime;
    }

    void Perf::perfSummerTime(){
        // std::cout <<"OUTPUT_DATA : " << OUTPUT_DATA << "  total = "<< total <<endl;
        if(!OUTPUT_DATA) return;
        double meanTime = totalTime / total;
        cout<<message;
        printf(" perf Timing total %d maxTiming:%f meanTiming:%f minTiming:%f\n", total, maxTime, meanTime, minTime);
    }