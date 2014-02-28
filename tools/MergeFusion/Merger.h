//
//  Merger.h
//  MergeFusion
//
//  Created by Zeng, Zheng/Sloan-Kettering Institute on 12/23/13.
//  Copyright (c) 2013 Zeng, Zheng/Sloan-Kettering Institute. All rights reserved.
//

#ifndef __MergeFusion__Merger__
#define __MergeFusion__Merger__


#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <string.h>
#include "common.h"
#include "StandardFusion.h"



using namespace std;

class Merger
{

public:
    Merger();

    ~Merger();

    void loadFusionFile(const char* filename, CALLER_TYPE caller_id, bool normalize_chrom, bool normalize_gene);
    
    void outputStandardFusion(const char* filename);

    void sortBy(string sorting_criteria);
    
    string locusToChr(string locus);
    
    void loadHugoWorker(string& chr_key, string& name_key, string &name_value);
    
    void loadHugoFile(const char* filename);
    
    void printHugoMap();
    
private:
    vector<StandardFusion> fusion_vec;
    map<pair<string, string>, map<string, int> > hugo_map;
};


#endif /* defined(__MergeFusion__Merger__) */
