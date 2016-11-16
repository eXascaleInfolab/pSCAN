#ifndef _UTILITY_H_
#define _UTILITY_H_

#ifdef _MSC_VER
	#define _CRT_SECURE_NO_WARNINGS
#endif

#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <queue>
#include <set>
// #include <unordered_set>

typedef unsigned  Id;
typedef int  Degree;  // ATTENTION: degrees are specified as int in the input file

#ifdef __unix__
#include <sys/time.h>
#endif

FILE *open_file(const char *file_name, const char *mode);

#endif
