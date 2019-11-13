//
//  common.cpp
//  MergeFusion
//
//  Created by Zeng, Zheng/Sloan-Kettering Institute on 12/30/13.
//  Copyright (c) 2013 Zeng, Zheng/Sloan-Kettering Institute. All rights reserved.
//

#include "common.h"


const char* CALLER_TYPE_STRING[] = {"CHIMERASCAN", "DEFUSE", "FUSIONCATCHER", "FUSIONHUNTER", "MAPSPLICE", "STAR"};
const bool CALLER_HEADER[] = {true, true, true, true, false, false};

void removeSpace(string& str)
{
    str.erase(remove_if(str.begin(), str.end(), ::isspace), str.end());
}


void split(const string& line, char delim, vector<std::string>& parsed_item, bool ignore_empty_item)
{
    stringstream my_ss(line);
    string item;
    while (getline(my_ss, item, delim))
    {
#ifdef _DEBUG
        cout << "[DEBUG]: Parsed item: " << item << endl;
#endif
        if (ignore_empty_item && item.empty())
            continue;// use to skip empty item
        parsed_item.push_back(item);
    }
}
