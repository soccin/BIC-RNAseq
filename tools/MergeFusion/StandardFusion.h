//
//  StandardFusion.h
//  MergeFusion
//
//  Created by Zeng, Zheng/Sloan-Kettering Institute on 12/23/13.
//  Copyright (c) 2013 Zeng, Zheng/Sloan-Kettering Institute. All rights reserved.
//

#ifndef __MergeFusion__StandardFusion__
#define __MergeFusion__StandardFusion_

#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <string.h>
#include "common.h"

using namespace std;


class StandardFusion
{
public:
    
    StandardFusion(string line, CALLER_TYPE _caller_id, bool bool_normalize_chrom, bool bool_normalize_gene,  map<pair<string, string>, map<string, int> >& hugo_map);
    
    ~StandardFusion();
    
    void init();
    
    void parseChimerascan(string line);
    
    void parseDefuse(string line);
    
    void parseFusioncatcher(string line);
    
    void parseFusionHunter(string line);
    
    void parseMapsplice(string line);
    
    void parseStar(string line);
    
    void normalize_chrom_name_worker(string& chrom_name);
    
    void normalize_chrom_name();
    
    void normalize_gene_name_worker(string& chrom_name, string &gene_name, map<pair<string, string>, map<string, int> >& hugo_map);
    
    void normalize_gene_name(map<pair<string, string>, map<string, int> >& hugo_map);
    
    bool operator<( const StandardFusion& other) const;
    
    friend ostream& operator<<(ostream &output, const StandardFusion& my_fusion);
    
    
    CALLER_TYPE caller_id;
    string chr1;
    string chr2;
    int start1;
    int start2;
    int end1;
    int end2;
    int breakpoint1;
    int breakpoint2;
    string strand1;
    string strand2;
    int total_coverage;
    int spanning_coverage;
    int breakpoint_coverage;
    string gene1;
    string gene2;
    string transcript1;
    string transcript2;
    string info;
};

void output_fusion_header(ostream &output);

bool sortByBreakpointCoverage(const StandardFusion &l, const StandardFusion &r);

bool sortByCaller(const StandardFusion &l, const StandardFusion &r);

bool sortByGenomicPos(const StandardFusion &l, const StandardFusion &r);

bool sortByGeneSymbol(const StandardFusion &l, const StandardFusion &r);








#endif /* defined(__MergeFusion__StandardFusion__) */
