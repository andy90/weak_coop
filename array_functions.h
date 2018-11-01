//
//  array_functions.h
//  WeakCoorp
//
//  Created by Ang Gao on 12/12/17.
//  Copyright © 2017 Ang Gao. All rights reserved.
//

#ifndef array_functions_h
#define array_functions_h
#include <iostream>
#include <vector>
#include <random>
#include <cstdio>
#include <cstdlib>
#include "names.h"
std::vector<double> array_normaldis(double mean, double sigma, int N, Rand_engine &generator);
std::vector<double> array_unidis(double minimum, double maximum, int N, Rand_engine &generator);
double norm(std::vector<double> &a, std::vector<double> &b);
double normsquare(std::vector<double> &a, std::vector<double> &b);
double array_sum(std::vector<double> &a);
std::vector<double> array_divide(std::vector<double> &a, double s);
void pbc(std::vector<double> &a, double box);
double normsquare2(std::vector<double> &a); 
double normcube(std::vector<double> &a, std::vector<double> &b);
#endif /* array_functions_h */
