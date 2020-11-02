#ifndef HAPCOLOR_H
#define HAPCOLOR_H

#include <iostream>
#include <list>
#include <vector>
#include <string>
#include <set>
#include <unordered_set>
#include <climits>
#include <iostream>
#include <memory>
#include <iterator>
#include <algorithm>
#include <sstream>
#include <tuple>
#include <fstream>
#include <math.h>
#include <cmath>
#include <map>
#include <cassert>
#include <time.h>
void compute_hapcolor(std::vector<int>*,std::vector<int>*, std::vector<int>*, std::vector<std::string>*, float, std::vector<int>*, std::vector<int>*, std::vector<float>*, std::string, std::string);
std::tuple<float, int, int> minimum(float*, int, bool, int, float);
float fwd_bkw(std::vector<int>, std::vector<std::vector<int>>, std::string);
float emission(int, int);
float transition(int, int);
float logSum(float, float);
#endif