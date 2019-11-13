//
//  common.h
//  MergeFusion
//
//  Created by Zeng, Zheng/Sloan-Kettering Institute on 12/30/13.
//  Copyright (c) 2013 Zeng, Zheng/Sloan-Kettering Institute. All rights reserved.
//

#ifndef __MergeFusion__common__
#define __MergeFusion__common__

#include <iostream>
#include <string.h>
#include <stdio.h>
#include <sstream>
#include <vector>
#include <map>
#include <stdlib.h>
#include <algorithm>

using namespace std;
//#define _DEBUG

enum CALLER_TYPE {CHIMERASCAN, DEFUSE, FUSIONCATCHER, FUSIONHUNTER, MAPSPLICE, STAR};
extern const char* CALLER_TYPE_STRING[]; //defined in StandardFusion.cpp
extern const bool CALLER_HEADER[]; //defined in StandardFusion.cpp

void removeSpace(string& str);

void split(const string& line, char delim, vector<std::string>& parsed_item, bool ignore_empty_item = false);

#endif /* defined(__MergeFusion__common__) */
