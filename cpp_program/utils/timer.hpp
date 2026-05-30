#ifndef timer_h
#define timer_h

#include <time.h>

struct Timer {
    clock_t start, stop;
    void start_timer() {start = clock();};
    void stop_timer() {stop = clock();};
    Timer() {start_timer();};
    float getTimeInSec(){return double(stop - start) / double(CLOCKS_PER_SEC);};
    float getTimeInMiliSec(){return double(stop - start) / double(CLOCKS_PER_SEC) * 1000.0;};
    float getTimeInMicroSec(){return double(stop - start) / double(CLOCKS_PER_SEC) * 1000.0 * 1000.0;};
    float getTimeInNanoSec(){return double(stop - start) / double(CLOCKS_PER_SEC) * 1000.0 * 1000.0 * 1000.0;};
};

#endif