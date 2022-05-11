#include <fenv.h>
#include "qcat_dataAnalysis.h"
#include <iostream>
#include "SZ3/api/sz.hpp"
#include "SZ3/utils/Config.hpp"
#include <stdlib.h>
#include <cmath>
#include <random>
#include <algorithm>
#include <mutex>
#include <thread>
#include <filesystem>
#include <vector>
#include <string>
#include <future>

#define OUTFILE "/home/ac.arhammkhan/sz3_analyze_data/sz3_data.csv"

namespace fs = std::filesystem;

int analyze_routine(const char* targetFile, size_t r1, size_t r2, size_t r3, size_t r4, size_t r5, size_t nbEle);
