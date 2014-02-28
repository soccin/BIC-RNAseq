//
//  main.cpp
//  MergeFusion
//
//  Created by Zeng, Zheng/Sloan-Kettering Institute on 12/23/13.
//  Copyright (c) 2013 Zeng, Zheng/Sloan-Kettering Institute. All rights reserved.
//

#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <string.h>
#include <getopt.h>
#include "common.h"
#include "Merger.h"




const string version = "MergeFusion 1.0.2";

vector< pair<string, CALLER_TYPE> > input_files;
string output_file;
string sorting_criteria = "GENOMIC";
bool normalize_chrom = true;
bool normalize_gene = false;
string hugo_data_file;



void printUsage(string msg = "")
{
    cout << endl;
    cout << version << endl;
    cout << "Usage: " << endl;
    cout << "\t[Each input file option can be specified multiple times]" << endl;
    cout << "\t--chimerascan        <string>                        Input chimerascan fusion file (chimeras.bedpe)" << endl;
    cout << "\t--defuse             <string>                        Input defuse fusion file (results.tsv)" << endl;
    cout << "\t--fusioncatcher      <string>                        Input fusioncatcher fusion file (final-list_candidate-fusion-genes.txt)" << endl;
    //cout << "\t--fusionhunter       <string>                        Input fusionhunter fusion file" << endl;
    cout << "\t--star               <string>                        Input star fusion file (Chimeric.out.junction)" << endl;
    cout << "\t--mapsplice          <string>                        Input mapsplice fusion file (fusions_well_annotated.txt)" << endl;
    cout << endl;
    cout << "\t--output             <string>                        Output fusion file to write the merged result" << endl;
    cout << "\t--sort               <string>                        Sort the result by either GENOMIC, CALLER, GENE, BREAKPOINT_COVERGAE, [Default=GENOMIC]" << endl;
    //cout << "\t--normalize_chrom    <int>                           Normalize the chromosome name to 'chr*'. 0=off, 1=on. [Default=1]" << endl;
    cout << "\t--normalize_gene     <string>                        Normalize the gene symbol to HUGO symbol with the hugo data file" << endl;
    cout << endl;
    
    if(!msg.empty())
        cerr << msg << endl;
    exit(1);
}

static struct option long_options[] =
{
    {"chimerascan",     required_argument,       0,     'C'},
    {"defuse",          required_argument,       0,     'D'},
    {"fusioncatcher",   required_argument,       0,     'F'},
    {"fusionhunter",    required_argument,       0,     'H'},
    {"star",            required_argument,       0,     'S'},
    {"mapsplice",       required_argument,       0,     'M'},
    {"output",          required_argument,       0,     'o'},
    {"sort",            required_argument,       0,     's'},
    //{"normalize_chrom", required_argument,       0,     'c'},
    {"normalize_gene",  required_argument,       0,     'g'},
    {0, 0, 0, 0}
};



void parseOption(int argc, char* argv[])
{
    if(argc == 1)
        printUsage();
    int next_option;
    int option_index = 0;
    do
    {
        next_option = getopt_long(argc, const_cast<char**>(argv), "C:D:F:H:S:M:o:s:g:", long_options, &option_index);
        switch(next_option)
        {
            case 'C':
                input_files.push_back(make_pair(optarg, CHIMERASCAN));
                break;
            case 'D':
                input_files.push_back(make_pair(optarg, DEFUSE));
                break;
            case 'F':
                input_files.push_back(make_pair(optarg, FUSIONCATCHER));
                break;
            case 'H':
                input_files.push_back(make_pair(optarg, FUSIONHUNTER));
                break;
            case 'S':
                input_files.push_back(make_pair(optarg, STAR));
                break;
            case 'M':
                input_files.push_back(make_pair(optarg, MAPSPLICE));
                break;
            case 'o':
                output_file = optarg;
                break;
            case 's':
                sorting_criteria = optarg;
                if(sorting_criteria != "CALLER" && sorting_criteria != "GENOMIC" && sorting_criteria != "GENE" && sorting_criteria != "BREAKPOINT_COVERAGE")
                    printUsage("[ERROR]: invalid value for --sort");
                break;
            /*case 'c':
                if(strcmp(optarg, "0")  == 0)
                    normalize_chrom = false;
                else if(strcmp(optarg, "1")  == 0)
                    normalize_chrom = true;
                else
                    printUsage("[ERROR]: invalid value for --normalize_chrom");
                break;*/
            case 'g':
                hugo_data_file = optarg;
                normalize_gene = true;
                break;
            case -1:
                break; //parsed all options
            default:
                printUsage(string("[ERROR]: Argument error: ") + argv[optind - 1]);
        }
    } while(next_option != -1);
    if(input_files.size() == 0)
        printUsage("[ERROR]: Please specify at least one input file");
    if(output_file.empty())
        printUsage("[ERROR]: Please specify output file");
    #ifdef _DEBUG
    cout << "[DEBUG]: Parsing options complete." << endl;
    #endif
}




int main(int argc, char* argv[])
{
    parseOption(argc, argv);
    Merger my_merger;
    if(normalize_gene)
        my_merger.loadHugoFile(hugo_data_file.c_str());
    //my_merger.printHugoMap();
    for(size_t i = 0; i < input_files.size(); i++)
    {
        my_merger.loadFusionFile(input_files[i].first.c_str(), input_files[i].second, normalize_chrom, normalize_gene);
    }
    if(!sorting_criteria.empty())
        my_merger.sortBy(sorting_criteria);
    my_merger.outputStandardFusion(output_file.c_str());
    return 0;
}

