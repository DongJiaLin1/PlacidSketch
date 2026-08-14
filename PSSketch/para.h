#ifndef PARA_H
#define PARA_H

#include <assert.h>
#include <time.h>


unordered_map<uint32_t, int> ground_truth_f;
unordered_map<uint32_t, int> ground_truth_p;
unordered_map<uint32_t, int> ground_truth_w;


//PISketch runtime parameters (set from MEM_CONFIG in the experiment loop)
// PI_size = 10.125 bytes per entry
extern int PI_X;
extern int PI_Y;
extern int BF_len;
extern int L;
extern int W_thr;
#define PI_size 10.125


//Common
#define P_thr 200
#define D_thr 3.0
#define WindowSize 1000
bool w_bit;

//Memory configurations (total sketch memory = PI_memory + BF_memory)
#define MEM_COUNT 5
#define MEM_400KB 0
#define MEM_500KB 1
#define MEM_600KB 2
#define MEM_700KB 3
#define MEM_800KB 4

//Memory configuration: {PI_X, PI_Y, BF_len}
// Target total ~400KB,500KB,600KB,700KB,800KB
const int MEM_CONFIG[MEM_COUNT][3] = {
    {394, 100, 10240},   // ~400KB  (PI=399KB + BF=10KB)
    {492, 100, 10240},   // ~500KB  (PI=498KB + BF=10KB)
    {591, 100, 10240},   // ~600KB  (PI=598KB + BF=10KB)
    {689, 100, 10240},   // ~700KB  (PI=698KB + BF=10KB)
    {788, 100, 10240},   // ~800KB  (PI=798KB + BF=10KB)
};
const char* MEM_NAMES[MEM_COUNT] = {
    "400KB", "500KB", "600KB", "700KB", "800KB"
};

//Others
#define Seed 99
#define TopN 1000
#define MAX_INSERT_PACKAGE 32000000

// Input path for the binary stream of flow keys (one uint32_t per record).
// Override at compile time, e.g.:
//   cmake -DDATA_BIN=/path/to/your/stream.bin ...
#ifndef DATA_BIN
#define DATA_BIN "stream.bin"
#endif

#define BF_size 1


//developer info
#define fmax 17057
#define pmax 2452


inline long pow(int a, int x){  //only for 2^x
    return 1L<<x;
}

struct report_item{
    uint32_t ID;
    int p;
    double D;
};



struct timespec pin_s, pin_e;
double Tsum=0;
int Tcnt=0;
inline void PIN(char p){
    if(p=='s')clock_gettime(CLOCK_MONOTONIC, &pin_s);
    else if(p=='e'){
        clock_gettime(CLOCK_MONOTONIC, &pin_e);
        long seconds = pin_e.tv_sec - pin_s.tv_sec;
        long nanoseconds = pin_e.tv_nsec - pin_s.tv_nsec;
        Tsum += (seconds + nanoseconds*1e-9);
        Tcnt ++;
    }
    else{  //'o'
        ;
    }
}




#endif