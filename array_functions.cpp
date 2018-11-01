//
//  array_functions.cpp
//  WeakCoorp
//
//  Created by Ang Gao on 12/12/17.
//  Copyright © 2017 Ang Gao. All rights reserved.
//

#include "array_functions.h"
#include <iostream>
#include <vector>
#include <random>
#include <cstdio>
#include <cstdlib>

std::vector<double> array_normaldis(double mean, double sigma, int N, Rand_engine &generator){
    std::vector<double> dis(N);
    std::normal_distribution<double> d(mean,sigma);
    for (int i=0; i<N; i++){
        dis[i]=d(generator);
    }
    
    return dis;
}

double norm(std::vector<double> &a, std::vector<double> &b){
    double d=0.0;
    if (a.size()==b.size()){
        for (int i=0; i<a.size(); i++){
            d+=pow(a[i]-b[i],2);
        }
    }else{
        std::cout<<"size mismatch\n";
        return -100;
    }
    d=sqrt(d);
    return d;
}

double normsquare(std::vector<double> &a, std::vector<double> &b){
    double d=0.0;
    if (a.size()==b.size()){
        for (int i=0; i<a.size(); i++){
            d+=(a[i]-b[i])*(a[i] - b[i]);
        }
    }else{
        std::cout<<"size mismatch\n";
        return -100;
    }
    return d;
}

double array_sum(std::vector<double> &a){
    double s=0.0;
    for (int i=0; i<a.size(); i++){
        s+=a[i];
    }
    return s;
}

std::vector<double> array_divide(std::vector<double> &a, double s){
    int n=a.size();
    std::vector<double> b(n);
    if (s!=0){
        for (int i=0; i<n; i++){
            b[i]=a[i]/s;
        }
    }
    return b;
}

void pbc(std::vector<double> &a, double box){
    int n=a.size();
    for (int i=0; i<n; i++){
        a[i]=a[i]-round(a[i]/(box/2))*box;
    }
}

std::vector<double> array_unidis(double minimum, double maximum, int N, Rand_engine &generator){
    std::vector<double> dis(N);
    std::uniform_real_distribution<double> uni(minimum, maximum);
    for (int i=0; i<N; i++) {
        dis[i]=uni(generator);
    }
    return dis;
}

double normsquare2(std::vector<double> &a) {
    double d = 0.0;

    for (int i = 0; i<a.size(); i++) {
        d += a[i] * a[i];
    }
    return d;
}

double normcube(std::vector<double> &a, std::vector<double> &b) {
    double d = 0.0;
    if (a.size() == b.size()) {
        for (int i = 0; i<a.size(); i++) {
            d += (a[i] - b[i])*(a[i] - b[i]);
        }
    }
    else {
        std::cout << "size mismatch\n";
        return -100;
    }
    d = d * sqrt(d);
    return d;
}

