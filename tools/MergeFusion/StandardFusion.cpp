//
//  StandardFusion.cpp
//  MergeFusion
//
//  Created by Zeng, Zheng/Sloan-Kettering Institute on 12/23/13.
//  Copyright (c) 2013 Zeng, Zheng/Sloan-Kettering Institute. All rights reserved.
//

#include "StandardFusion.h"

StandardFusion::StandardFusion(string line, CALLER_TYPE _caller_id, bool bool_normalize_chrom, bool bool_normalize_gene, map<pair<string, string>, map<string, int> >& hugo_map)
{
    init();
    caller_id = _caller_id;
    switch(caller_id)
    {
        case CHIMERASCAN: parseChimerascan(line); break;
        case DEFUSE: parseDefuse(line); break;
        case FUSIONCATCHER: parseFusioncatcher(line); break;
        case FUSIONHUNTER: parseFusionHunter(line); break;
        case MAPSPLICE: parseMapsplice(line); break;
        case STAR: parseStar(line); break;
        default: cerr << "Unknown fusion file type: " << caller_id << endl; exit(1);
    }
    if(bool_normalize_chrom)
        normalize_chrom_name();
    if(bool_normalize_gene)
        normalize_gene_name(hugo_map);
}

StandardFusion::~StandardFusion() {}

void StandardFusion::init()
{
    start1 = -1;
    start2 = -1;
    end1 = -1;
    end2 = -1;
    breakpoint1 = -1;
    breakpoint2 = -1;
    total_coverage = -1;
    spanning_coverage = -1;
    breakpoint_coverage = -1;
}


void StandardFusion::parseChimerascan(string line) //done
{
    vector<string> parsed_item;
    vector<string> breakpoint_reads;
    split(line, '\t', parsed_item);
    chr1 = parsed_item[0];
    start1 = atoi(parsed_item[1].c_str()) + 1;
    end1 = atoi(parsed_item[2].c_str()) + 1;
    chr2 = parsed_item[3];
    start2 = atoi(parsed_item[4].c_str()) + 1;
    end2 = atoi(parsed_item[5].c_str()) + 1;
    strand1 = parsed_item[8];
    strand2 = parsed_item[9];
    breakpoint1 = (strand1 == "+" ? end1 : start1);
    breakpoint2 = (strand2 == "+" ? start2 : end2);
    gene1 = parsed_item[12];
    gene2 = parsed_item[13];
    transcript1 = parsed_item[10];
    transcript2 = parsed_item[11];
    total_coverage = atoi(parsed_item[16].c_str());
    spanning_coverage = atoi(parsed_item[17].c_str());
    split(parsed_item[21], '>', breakpoint_reads, true);
    breakpoint_coverage = static_cast<int>(breakpoint_reads.size());
    caller_id = CHIMERASCAN;
    
    info.append("CHIMERA_CLUSTER_ID=");info.append(parsed_item[6]);
    info.append(";SCORE=");info.append(parsed_item[7]);
    info.append(";TYPE=");info.append(parsed_item[14]);
    info.append(";DISTANCE=");info.append(parsed_item[15]);
    info.append(";UNIQUE_ALIGNMENT_POSITIONS=");info.append(parsed_item[18]);
    info.append(";ISOFORM_FRACTION_5P=");info.append(parsed_item[19]);
    info.append(";ISOFORM_FRACTION_3P=");info.append(parsed_item[20]);
    info.append(";CHIMERA_IDS="); info.append(parsed_item[22]);
}

void StandardFusion::parseDefuse(string line) // problem with strand, not consistent with other callers
{
    vector<string> parsed_item;
    split(line, '\t', parsed_item);
    chr1 = parsed_item[24];
    chr2 = parsed_item[25];
    start1 = atoi(parsed_item[32].c_str());
    start2 = atoi(parsed_item[33].c_str());
    end1 = atoi(parsed_item[26].c_str());
    end2 = atoi(parsed_item[27].c_str());
    breakpoint1 = atoi(parsed_item[37].c_str());
    breakpoint2 = atoi(parsed_item[38].c_str());
    strand1 = parsed_item[39];
    strand2 = parsed_item[40];
    //strand1 = parsed_item[22];
    //strand2 = parsed_item[23];
    gene1 = parsed_item[30];
    gene2 = parsed_item[31];
    transcript1 = "NA";
    transcript2 = "NA";
    spanning_coverage = atoi(parsed_item[56].c_str());
    breakpoint_coverage = atoi(parsed_item[2].c_str());
    total_coverage = spanning_coverage + breakpoint_coverage;
    caller_id = DEFUSE;
    
    info.append("CLUSTER_ID=");info.append(parsed_item[0]);
    info.append(";SPLITR_SPAN_PVALUE="); info.append(parsed_item[3]);
    info.append(";SPLITR_POS_PVALUE="); info.append(parsed_item[4]);
    info.append(";SPLITR_MIN_PVALUE="); info.append(parsed_item[5]);
    info.append(";ADJACENT="); info.append(parsed_item[6]);
    info.append(";ALTSPLICE="); info.append(parsed_item[7]);
    info.append(";BREAK_ADJ_ENTROPY1="); info.append(parsed_item[8]);
    info.append(";BREAK_ADJ_ENTROPY2="); info.append(parsed_item[9]);
    info.append(";BREAK_ADJ_ENTROPY_MIN="); info.append(parsed_item[10]);
    info.append(";BREAKPOINT_HOMOLOGY="); info.append(parsed_item[11]);
    info.append(";BREAKSEQS_ESTISLANDS_PERCIDENT="); info.append(parsed_item[12]);
    info.append(";CDNA_BREAKSEQS_PERCIDENT="); info.append(parsed_item[13]);
    info.append(";DELETION="); info.append(parsed_item[14]);
    info.append(";EST_BREAKSEQS_PERCIDENT="); info.append(parsed_item[15]);
    info.append(";EVERSION="); info.append(parsed_item[16]);
    info.append(";EXONBOUNDARIES="); info.append(parsed_item[17]);
    info.append(";EXPRESSION1="); info.append(parsed_item[18]);
    info.append(";EXPRESS2="); info.append(parsed_item[19]);
    info.append(";GENE_LOCATION1="); info.append(parsed_item[28]);
    info.append(";GENE_LOCATION2="); info.append(parsed_item[29]);
    info.append(";GENOME_BREAKSEQS_PERCIDENT="); info.append(parsed_item[36]);
    info.append(";INTERCHROMOSOMAL="); info.append(parsed_item[41]);
    info.append(";INTERRUPTED_INDEX1="); info.append(parsed_item[42]);
    info.append(";INTERRUPTED_INDEX2="); info.append(parsed_item[43]);
    info.append(";INVERSION="); info.append(parsed_item[44]);
    info.append(";LIBRARY_NAME="); info.append(parsed_item[45]);
    info.append(";MAX_MAP_COUNT="); info.append(parsed_item[46]);
    info.append(";MAX_REPEAT_PROPORTION="); info.append(parsed_item[47]);
    info.append(";MEAN_MAP_COUNT="); info.append(parsed_item[48]);
    info.append(";MIN_MAP_COUNT="); info.append(parsed_item[49]);
    info.append(";NUM_MULTI_MAP="); info.append(parsed_item[50]);
    info.append(";NUM_SPLICE_VARIANTS="); info.append(parsed_item[51]);
    info.append(";ORF="); info.append(parsed_item[52]);
    info.append(";READ_THROUGH="); info.append(parsed_item[53]);
    info.append(";REPEAT_PROPORTION1="); info.append(parsed_item[54]);
    info.append(";REPEAT_PROPORTION2="); info.append(parsed_item[55]);
    info.append(";SPAN_COVERAGE1="); info.append(parsed_item[57]);
    info.append(";SPAN_COVERAGE2="); info.append(parsed_item[58]);
    info.append(";SPAN_COVERAGE_MAX="); info.append(parsed_item[59]);
    info.append(";SPAN_COVERAGE_MIN"); info.append(parsed_item[60]);
    info.append(";SPLICE_SCORE="); info.append(parsed_item[61]);
    info.append(";SPLICING_INDEX1="); info.append(parsed_item[62]);
    info.append(";SPLICING_INDEX2="); info.append(parsed_item[63]);
}

void StandardFusion::parseFusioncatcher(string line) //done
{
    vector<string> parsed_item;
    vector<string> pos_parsed_item1;
    vector<string> pos_parsed_item2;
    split(line, '\t', parsed_item);
    split(parsed_item[8], ':', pos_parsed_item1);
    split(parsed_item[9], ':', pos_parsed_item2);
    chr1 = pos_parsed_item1[0];
    chr2 = pos_parsed_item2[0];
    breakpoint1 = atoi(pos_parsed_item1[1].c_str());
    breakpoint2 = atoi(pos_parsed_item2[1].c_str());
    strand1 = pos_parsed_item1[2];
    strand2 = pos_parsed_item2[2];
    gene1 = parsed_item[0];
    gene2 = parsed_item[1];
    transcript1 = "NA";
    transcript2 = "NA";
    start1 = -1;
    start2 = -1;
    end1 = -1;
    end2 = -1;
    spanning_coverage = atoi(parsed_item[4].c_str());
    breakpoint_coverage = atoi(parsed_item[5].c_str());
    total_coverage = spanning_coverage + breakpoint_coverage;
    caller_id = FUSIONCATCHER;
    
    info.append("COUNT_COMMON_MAPPING_READS="); info.append(parsed_item[3]);
    info.append(";LONGEST_ANCHOR="); info.append(parsed_item[6]);
    info.append(";FUSION_FINDING_METHOD="); info.append(parsed_item[7]);
    info.append(";GENE1_ID="); info.append(parsed_item[10]);
    info.append(";GENE2_ID="); info.append(parsed_item[11]);
    info.append(";EXON1_ID="); info.append(parsed_item[12]);
    info.append(";EXON2_ID="); info.append(parsed_item[13]);
}

void StandardFusion::parseFusionHunter(string line)
{
    caller_id = FUSIONHUNTER;
}

void StandardFusion::parseMapsplice(string line) // done
{
    vector<string> parsed_item;
    vector<string> chr_parsed_item;
    split(line, '\t', parsed_item);
    split(parsed_item[0], '~', chr_parsed_item);
    chr1 = chr_parsed_item[0];
    chr2 = chr_parsed_item[1];
    strand1 = parsed_item[5][0];
    strand2 = parsed_item[5][1];
    breakpoint1 = atoi(parsed_item[1].c_str());
    breakpoint2 = atoi(parsed_item[2].c_str());
    start1 = (strand1 == "+" ? atoi(parsed_item[28].c_str()) : breakpoint1);
    end1 = (strand1 == "-" ? atoi(parsed_item[28].c_str()) : breakpoint1);
    start2 = (strand2 == "-" ? atoi(parsed_item[29].c_str()) : breakpoint2);
    end2 = (strand2 == "+" ? atoi(parsed_item[29].c_str()) : breakpoint2);
    gene1 = parsed_item[60];
    gene1 = gene1.substr(0, gene1.length() - 1);
    gene2 = parsed_item[61];
    gene2 = gene2.substr(0, gene2.length() - 1);
    transcript1 = "NA";
    transcript2 = "NA";
    breakpoint_coverage = atoi(parsed_item[4].c_str());
    spanning_coverage = atoi(parsed_item[27].c_str());
    total_coverage = spanning_coverage + breakpoint_coverage;
    caller_id = MAPSPLICE;
    
    info.append("ID="); info.append(parsed_item[3]);
    info.append(";ENTROPY="); info.append(parsed_item[10]);
    info.append(";FLANK_STRING="); info.append(parsed_item[12]);
    info.append(";MIN_MISMATCH="); info.append(parsed_item[13]);
    info.append(";MAX_MISMATCH="); info.append(parsed_item[14]);
    info.append(";AVE_MISMATCH="); info.append(parsed_item[15]);
    info.append(";MAX_MIN_SUFFIX="); info.append(parsed_item[16]);
    info.append(";MAX_MIN_PREFIX="); info.append(parsed_item[17]);
    info.append(";MIN_ANCHOR_DIFF="); info.append(parsed_item[18]);
    info.append(";UNIQUE_BREAKPOINT_COVERAGE="); info.append(parsed_item[19]);
    info.append(";MULTI_BREAKPOINT_COVERAGE="); info.append(parsed_item[20]);
    info.append(";PAIRED_READ_COUNT="); info.append(parsed_item[21]);
    info.append(";SINGLE_READ_COUNT="); info.append(parsed_item[26]);
    info.append(";DONOR_ISOFORMS="); info.append(parsed_item[30]);
    info.append(";ACCEPTOR_ISOFORMS="); info.append(parsed_item[31]);
    info.append(";MIN_DONOR_ISOFORM_LENGTH="); info.append(parsed_item[36]);
    info.append(";MAX_DONOR_ISOFORM_LENGTH="); info.append(parsed_item[37]);
    info.append(";MIN_ACCEPTOR_ISOFORM_LENGTH="); info.append(parsed_item[38]);
    info.append(";MAX_ACCEPTOR_ISOFORM_LENGTH="); info.append(parsed_item[39]);
    info.append(";DONOR_MATCH_NORMAL="); info.append(parsed_item[52]);
    info.append(";ACCEPTOR_MATCH_NORMAL="); info.append(parsed_item[53]);
    info.append(";MATCH_GENE_STRAND="); info.append(parsed_item[56]);
    info.append(";ANNOTATE_TYPE="); info.append(parsed_item[57]);
    info.append(";FUSION_TYPE="); info.append(parsed_item[58]);
}

void StandardFusion::parseStar(string line)  //done
{
    vector<string> parsed_item;
    split(line, '\t', parsed_item);
    chr1 = parsed_item[0];
    chr2 = parsed_item[3];
    start1 = -1;
    start2 = -1;
    end1 = -1;
    end2 = -1;
    breakpoint1 = atoi(parsed_item[1].c_str());
    breakpoint2 = atoi(parsed_item[4].c_str());
    strand1 = parsed_item[2];
    strand2 = parsed_item[5];
    breakpoint_coverage = 1;
    spanning_coverage = -1;
    total_coverage = -1;
    gene1 = "NA";
    gene2 = "NA";
    transcript1 = "NA";
    transcript2 = "NA";
    caller_id = STAR;
    
    info.append("JUNCTION_TYPE="); info.append(parsed_item[6]);
    info.append(";REPEAT_LEFT="); info.append(parsed_item[7]);
    info.append(";REPEAT_RIGHT="); info.append(parsed_item[8]);
}

void StandardFusion::normalize_chrom_name_worker(string& chrom_name)
{
    if(chrom_name == "M" || chrom_name == "Mt" || chrom_name == "chrM" || chrom_name == "chrMt")
        chrom_name = "MT";
    if(chrom_name.substr(0, 3) != "chr")
        chrom_name = "chr" + chrom_name;
}

void StandardFusion::normalize_chrom_name()
{
    normalize_chrom_name_worker(chr1);
    normalize_chrom_name_worker(chr2);
}

void StandardFusion::normalize_gene_name_worker(string& chrom_name, string& gene_name, map<pair<string, string>, map<string, int> >& hugo_map)
{
    map<string, int> hugo_names;
    vector<string> parsed_item;
    split(gene_name, ',', parsed_item, true);
    for(size_t i = 0 ; i < parsed_item.size(); i++)
    {
        string ori_sub_gene_name = parsed_item[i];
        string sub_gene_name = parsed_item[i];
        removeSpace(sub_gene_name);
        transform(sub_gene_name.begin(), sub_gene_name.end(),sub_gene_name.begin(), ::toupper);
        pair<string, string> combined_key = make_pair(chrom_name, sub_gene_name);
        string new_hugo_name;
        if(hugo_map.find(combined_key) != hugo_map.end() && hugo_map[combined_key].size() == 1)
        {
            new_hugo_name = hugo_map[combined_key].begin()->first;
        }
        else
        {
            new_hugo_name = ori_sub_gene_name;
        }
        if(hugo_names.find(new_hugo_name) == hugo_names.end())
            hugo_names.insert(make_pair(new_hugo_name, 1));
        else
            hugo_names[new_hugo_name]++;
    }
    string final_name;
    for(map<string, int>::iterator it = hugo_names.begin(); it != hugo_names.end(); it ++)
    {
        if(!final_name.empty())
            final_name.append(",");
        final_name.append(it->first);
    }
    gene_name = final_name;
}

void StandardFusion::normalize_gene_name(map<pair<string, string>, map<string, int> >& hugo_map)
{
    if(gene1 != "NA")
        normalize_gene_name_worker(chr1, gene1, hugo_map);
    if(gene2 != "NA")
        normalize_gene_name_worker(chr2, gene2, hugo_map);
}

bool StandardFusion::operator<(const StandardFusion& other) const
{
    if(chr1 != other.chr1)
        return chr1 < other.chr1;
    if(strand1 != other.strand1)
        return strand1 < other.strand1;
    if(breakpoint1 != other.breakpoint1)
        return breakpoint1 < other.breakpoint1;
    if(chr2 != other.chr2)
        return chr2 < other.chr2;
    if(strand2 != other.strand2)
        return strand2 < other.strand2;
    if(breakpoint2 != other.breakpoint2)
        return breakpoint2 < other.breakpoint2;
    else
        return caller_id < other.caller_id;
}

void output_fusion_header(ostream &output)
{
    output << "CALLER_ID" << "\t" \
    << "GENE1"<< "\t" \
    << "GENE2" << "\t" \
    << "CHR1" << "\t" \
    << "STRAND1" << "\t" \
    << "BREAK_POINT1" << "\t" \
    << "CHR2" << "\t" \
    << "STRAND2" << "\t" \
    << "BREAK_POINT2" << "\t" \
    << "BREAKPOINT_COVERAGE" << "\t" \
    << "SPANNING_COVERAGE" << "\t" \
    << "TOTAL_COVERAGE" << "\t" \
    << "START1" << "\t" \
    << "END1" << "\t" \
    << "START2" << "\t" \
    << "END2" << "\t" \
    << "TRANSCRIPT1" << "\t" \
    << "TRANSCRIPT2" << "\t" \
    << "INFO" \
    << endl;
}


ostream& operator<<(ostream &output, const StandardFusion& my_fusion)
{
    output << CALLER_TYPE_STRING[my_fusion.caller_id] << "\t" \
    << my_fusion.gene1 << "\t" \
    << my_fusion.gene2 << "\t" \
    << my_fusion.chr1 << "\t" \
    << my_fusion.strand1 << "\t" \
    << my_fusion.breakpoint1 << "\t" \
    << my_fusion.chr2 << "\t" \
    << my_fusion.strand2 << "\t" \
    << my_fusion.breakpoint2 << "\t";
    if(my_fusion.breakpoint_coverage >= 0)
        output << my_fusion.breakpoint_coverage << "\t";
    else
        output << "NA" << "\t";
    if(my_fusion.spanning_coverage >= 0)
        output << my_fusion.spanning_coverage << "\t";
    else
        output << "NA" << "\t";
    if(my_fusion.total_coverage >= 0)
        output << my_fusion.total_coverage << "\t";
    else
        output << "NA" << "\t";
    if(my_fusion.start1 >= 0)
        output << my_fusion.start1 << "\t";
    else
        output << "NA" << "\t";
    if(my_fusion.end1 >= 0)
        output << my_fusion.end1 << "\t";
    else
        output << "NA" << "\t";
    if(my_fusion.start2 >= 0)
        output << my_fusion.start2 << "\t";
    else
        output << "NA" << "\t";
    if(my_fusion.end2 >= 0)
        output << my_fusion.end2 << "\t";
    else
        output << "NA" << "\t";
    output << my_fusion.transcript1 << "\t" \
    << my_fusion.transcript2 << "\t" \
    << my_fusion.info;
    return output;
}


bool sortByBreakpointCoverage(const StandardFusion &l, const StandardFusion &r)
{
    return l.breakpoint_coverage > r.breakpoint_coverage;
}

bool sortByCaller(const StandardFusion &l, const StandardFusion &r)
{
    return l.caller_id < r.caller_id;
}

bool sortByGenomicPos(const StandardFusion &l, const StandardFusion &r)
{
    if(l.chr1 != r.chr1)
        return l.chr1 < r.chr1;
    if(l.chr2 != r.chr2)
        return l.chr2 < r.chr2;
    if(l.breakpoint1 != r.breakpoint1)
        return l.breakpoint1 < r.breakpoint1;
    if(l.breakpoint2 != r.breakpoint2)
        return l.breakpoint2 < r.breakpoint2;
    if(l.strand1 != r.strand1)
        return l.strand1 < r.strand1;
    else
        return l.strand2 < r.strand2;
}

bool sortByGeneSymbol(const StandardFusion &l, const StandardFusion &r)
{
    if(l.gene1 != r.gene1)
        return l.gene1 < r.gene1;
    else
        return l.gene2 < r.gene2;
}



















