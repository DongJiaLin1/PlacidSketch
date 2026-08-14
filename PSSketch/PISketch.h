
#include"BOBHash32.h"
#include<queue>
#include<unordered_map>
#include<iostream>
#include<map>
#include<ctime>
#include<cstring>
#include<stdlib.h>
#include<vector>
using namespace std;

class BloomFilter{
    int k,m;
    BOBHash32 *bh1,*bh2,*bh3;
    bool *S;
public:
    BloomFilter(int _k,int _m):k(_k),m(_m){
        bh1 = new BOBHash32(11); 
        bh2 = new BOBHash32(17);
		bh3 = new BOBHash32(23);
        S = new bool[m];
        memset(S,0,sizeof(bool)*m);
    }
    ~BloomFilter(){
        delete []S;
    }
    void insert(uint32_t x){
    	int pos1 = bh1->run((const char*)&x,4) % m;
		int pos2 = bh2->run((const char*)&x,4) % m;
		int pos3 = bh3->run((const char*)&x,4) % m;
        S[pos1]=1;
		S[pos2]=1;
		S[pos3]=1;
    }
    bool query(uint32_t x){
    	int pos1 = bh1->run((const char*)&x,4) % m;
		int pos2 = bh2->run((const char*)&x,4) % m;
		int pos3 = bh3->run((const char*)&x,4) % m;
		return S[pos1]&&S[pos2]&&S[pos3];
    }
    void clear(){
        memset(S,0,sizeof(bool)*m);
    }
    int getMem(){
        return 1*m;
    }
};

struct Bucket{
    int *w;
	uint16_t *p;
    uint32_t *fp;
    bool *valid;
    int maxnum;
    Bucket(){}
    void init(int maxnum_){
        maxnum = maxnum_;
        valid = new bool[maxnum];
        memset(valid, 0, maxnum * sizeof(bool));
        w = new int[maxnum];
        p = new uint16_t[maxnum];
        fp = new uint32_t[maxnum];
    }
    void insert(uint32_t d,int wi){  
        int mi = -1;
        for(int i=0;i<maxnum;++i){
            if(!valid[i]) continue;
            if(fp[i] == d){ 
                w[i] += wi;    
                if(wi > 0) p[i] ++;  
                return;
            }
            if(mi < 0 || w[i] < w[mi]) mi = i;  
        }
        for(int i=0;i<maxnum;++i){ 
            if(!valid[i]){ 
                fp[i] = d;
                w[i] = wi;  
                p[i] = 1;  
                valid[i] = 1;
                return;
            }
        }
        w[mi] --;  
       
        if(w[mi] < 0){  
            fp[mi] = d; 
            w[mi] = wi + 1;
            p[mi] = 1;
        }
    }

	
};
class PISketch{
public:
    Bucket *bucket;
    int bucket_num;
    int num_per_bucket;
    BOBHash32 *bh;
    PISketch(int bn,int num_per_bucket_):bucket_num(bn){
        num_per_bucket = num_per_bucket_;
        bucket = new Bucket[bucket_num];
        for(int i=0;i<bucket_num;++i) bucket[i].init(num_per_bucket);
        bh = new BOBHash32(37);
    }
    void insert(uint32_t d,int wi){
        int pos = bh->run((const char*)&d,4) %bucket_num;
        bucket[pos].insert(d,wi);
    }


	uint32_t getfp(int x,int y){
		return bucket[x].fp[y];		
	}

	int getw(int x,int y){
		return bucket[x].w[y];		
	}

};