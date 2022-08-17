#pragma once

#include <cmath>
#include <cstdio>
#include <vector>
#include <iostream>
#include <algorithm>
#include <iomanip>

double pearsoncoeff(std::vector<uint> X, std::vector<uint> Y);
uint sum(std::vector<uint> a);
double mean(std::vector<uint> a);
uint sqsum(std::vector<uint> a);
double stdev(std::vector<uint> nums);
std::vector<uint> operator-(std::vector<uint> a, uint b);
std::vector<uint> operator*(std::vector<uint> a, std::vector<uint> b);
