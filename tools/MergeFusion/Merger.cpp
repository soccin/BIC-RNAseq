//
//  Merger.cpp
//  MergeFusion
//
//  Created by Zeng, Zheng/Sloan-Kettering Institute on 12/23/13.
//  Copyright (c) 2013 Zeng, Zheng/Sloan-Kettering Institute. All rights reserved.
//

#include "Merger.h"


Merger::Merger()
{}

Merger::~Merger()
{}


void Merger::loadFusionFile(const char* filename, CALLER_TYPE caller_id, bool normalize_chrom, bool normalize_gene)
{
#ifdef _DEBUG
    cout << "[DEBUG]: Loading " << CALLER_TYPE_STRING[caller_id] << " file: " << filename << endl;
#endif
    
    ifstream input_fs(filename);
    if(!input_fs)
    {
        cerr << "Error: fail to open input file: " << filename << endl;
        exit(1);
    }
    map<StandardFusion, size_t> collapse_star;
    string line;
    if(CALLER_HEADER[caller_id])
    {
        getline(input_fs, line);
#ifdef _DEBUG
        cout << "[DEBUG]: Header line: " << line << endl;
#endif
    }
    while(getline(input_fs, line))
    {
#ifdef _DEBUG
        cout << "[DEBUG]: Data line: " << line << endl;
#endif
        StandardFusion new_fusion_entry(line, caller_id, normalize_chrom, normalize_gene, hugo_map);
        if(caller_id != STAR)   // caller other than star
            fusion_vec.push_back(new_fusion_entry);
        else                    // collapse star entries for each fusion junction
        {
            if(collapse_star.find(new_fusion_entry) == collapse_star.end())
            {
                fusion_vec.push_back(new_fusion_entry);
                collapse_star.insert(make_pair(new_fusion_entry, fusion_vec.size() - 1));
            }
            else
            {
                size_t entry_index = collapse_star[new_fusion_entry];
                fusion_vec[entry_index].breakpoint_coverage ++;
            }
        }
    }
    input_fs.close();
}


void Merger::outputStandardFusion(const char* filename)
{
#ifdef _DEBUG
    cout << "[DEBUG]: Printing result to: " << filename << endl;
#endif
    ofstream output_fs(filename);
    if(!output_fs)
    {
        cerr << "Error: fail to open output file: " << filename << endl;
        exit(1);
    }
    output_fusion_header(output_fs);
    for(size_t i = 0; i < fusion_vec.size(); i++)
        output_fs << fusion_vec[i] << endl;
    output_fs.close();
}


void Merger::sortBy(string sorting_criteria)
{
    if(sorting_criteria == "CALLER")
    {
        sort(fusion_vec.begin(), fusion_vec.end(), sortByCaller);
    }
    else if(sorting_criteria == "GENOMIC")
    {
        sort(fusion_vec.begin(), fusion_vec.end(), sortByGenomicPos);
    }
    else if(sorting_criteria == "GENE")
    {
        sort(fusion_vec.begin(), fusion_vec.end(), sortByGeneSymbol);
    }
    else if(sorting_criteria == "BREAKPOINT_COVERAGE")
    {
        sort(fusion_vec.begin(), fusion_vec.end(), sortByBreakpointCoverage);
    }
    else
    {
        cerr << "[ERROR]: unrecognized sorting criteria '" << sorting_criteria << "', output won't be sorted" << endl;
        
    }
}


string Merger::locusToChr(string locus)
{
#ifdef _DEBUG
    cerr << "[DEBUG]: locus= " << locus << endl;
#endif
    if(locus.empty() || locus == "reserved")
        return "chrNA";
    if(locus == "mitochondria")
        return "chrMT";
    else
    {
        vector<string> p_split_item;
        split(locus, 'p', p_split_item);
        vector<string> q_split_item;
        split(p_split_item[0], 'q', q_split_item);
        vector<string> c_split_item;
        split(q_split_item[0], 'c', c_split_item);
        return "chr" + c_split_item[0];
    }
}

void Merger::loadHugoWorker(string& chr_key, string& name_key, string &name_value)
{
    vector<string> parsed_item;
    split(name_key, ',', parsed_item, true);
    for(size_t i = 0; i < parsed_item.size(); i++)
    {
        string sub_name_key = parsed_item[i];
        removeSpace(sub_name_key);
        transform(sub_name_key.begin(), sub_name_key.end(),sub_name_key.begin(), ::toupper);
        pair<string, string> combined_key= make_pair(chr_key, sub_name_key);
        if(hugo_map.find(combined_key) == hugo_map.end())
        {
            map<string, int> new_val_map;
            hugo_map.insert(make_pair(combined_key, new_val_map));
        }
        map<string, int>& val_map = hugo_map[combined_key];
        if(val_map.find(name_value) == val_map.end())
            val_map.insert(make_pair(name_value, 1));
        else
            val_map[name_value] ++;
    }
}

void Merger::loadHugoFile(const char* filename)
{
    int approved_symbol_index = 1;
    int approved_status_index = 3;
    int previous_symbol_index = 5;
    int synonyms_index = 7;
    int locus_index = 9;
    int accesion_index = 14;
    int refseq_id_index = 21;
    int vega_id_index = 25;
    
    // int entrez_id_ncbi_index = 28;
    // in omim_id_ncbi_index = 29;
    int refseq_id_ncbi_index = 30;
    //int uniprot_id_ncbi_index = 31;
    int ensembl_id_index = 32;
    int uscs_id_index = 33;
    ifstream input_fs(filename);
    if(!input_fs)
    {
        cerr << "Error: fail to open Hugo file: " << filename << endl;
        exit(1);
    }
    string line;
    getline(input_fs, line);
    while(getline(input_fs, line))
    {
        vector<string> parsed_item;
        split(line, '\t', parsed_item);
        if(parsed_item[approved_status_index] != "Approved")
            continue;
        string chr_key = locusToChr(parsed_item[locus_index]);
        if(parsed_item.size() >= 2)
            loadHugoWorker(chr_key, parsed_item[approved_symbol_index], parsed_item[approved_symbol_index]);
        if(parsed_item.size() >= 6)
            loadHugoWorker(chr_key, parsed_item[previous_symbol_index], parsed_item[approved_symbol_index]);
        if(parsed_item.size() >= 8)
            loadHugoWorker(chr_key, parsed_item[synonyms_index], parsed_item[approved_symbol_index]);
        if(parsed_item.size() >= 15)
            loadHugoWorker(chr_key, parsed_item[accesion_index], parsed_item[approved_symbol_index]);
        if(parsed_item.size() >= 22)
            loadHugoWorker(chr_key, parsed_item[refseq_id_index], parsed_item[approved_symbol_index]);
        if(parsed_item.size() >= 26)
            loadHugoWorker(chr_key, parsed_item[vega_id_index], parsed_item[approved_symbol_index]);
        if(parsed_item.size() >= 31)
            loadHugoWorker(chr_key, parsed_item[refseq_id_ncbi_index], parsed_item[approved_symbol_index]);
        if(parsed_item.size() >= 33)
            loadHugoWorker(chr_key, parsed_item[ensembl_id_index], parsed_item[approved_symbol_index]);
        if(parsed_item.size() >= 34)
            loadHugoWorker(chr_key, parsed_item[uscs_id_index], parsed_item[approved_symbol_index]);
    }
    input_fs.close();
}

void Merger::printHugoMap()
{
    map<pair<string, string>, map<string, int> >::iterator it;
    for(it = hugo_map.begin(); it != hugo_map.end(); it ++)
    {
        cout << it->first.first << "\t" << it->first.second << "\t" << it->second.size() << endl;
        /*if(it->second.size() > 1)
        {
            map<string, int>::iterator it2;
            for(it2 = it->second.begin(); it2 != it->second.end(); it2 ++)
                cout << it2->first << "\t" << it2->second << endl;
        }*/
    }
}



























