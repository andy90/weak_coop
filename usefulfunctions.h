#pragma once
#include <iostream>
#include <vector>
#include <map>
#include <list>
#include <ctime>
#include <algorithm>
#include <fstream>
#include <string>
#include <random>
#include "names.h"
#include "Organism.h"

double energy(std::vector<double> comp, std::vector<double> task);

void calc_ntE_ntE2(CompCluster_Energy &CCE, std::vector<double> &ntE, std::vector<double> &ntE2);

double costfunction(int N, std::vector<double> ntE, std::vector<double> ntE2);

void add_comp(Comp_Task &components, CompCluster_Energy &CCE, std::vector<double> new_comp, int new_ic, std::map<int, double> energy_ic, coopflag cflag);

void delete_comp(CompCluster_Energy &CCE, Comp_info &tasks_for_comp, Comp_Task &components, int ic, coopflag cflag);

void generate_mutatelist(Rand_engine &generator, std::vector<int> &mlist);

std::map<int, double> calc_tasksenergys(std::vector<double> &component, Comp_Task &tasks);

void mutate(std::list<Organism> &organisms, Rand_engine &generator);

void resample(std::list<Organism> &organisms, Rand_engine &generator);

void mutateresample(std::list<Organism> &organisms, Rand_engine &generator, std::ofstream &f); // do mutation and resample in one storke, and also write out the death rate of the mutated organism after selection

void write_completion_status(std::list<Organism> &organisms, std::ofstream &f);

void write_edge_percent(CompCluster_Energy &CCE, Comp_info &tasks_for_comp, std::vector<double> &ntE, std::vector<double> &ntE2, std::vector<double> &n_task, int norg);

bool nonzero(double a);

void write_mutations(std::list<Organism> &organisms, std::ofstream &f);

std::vector<double> write_ntE_av(std::list<Organism> &organisms, int M, std::ofstream &f);

void write_traj_mostprob(std::list<Organism> &organisms, std::ofstream &f, int isteps);

void write_traj(Organism &the_organism, std::ofstream &traj, int isteps);

void writeN(std::list<Organism> &organisms, std::ofstream &f);

void write_cost_av(std::list<Organism> &organisms, std::ofstream &f2);

std::vector<double> generate_delta(double d, int D, Rand_engine &generator_fixed);

void write_cluster_size(CompCluster_Energy &CCE, std::vector<double> &Ndis, int norg);

void write_csize_av(std::list<Organism> &organisms, std::ofstream &f);

void write_cross_dis(CompCluster_Energy &CCE, int norg, Comp_Task &tasks, std::ofstream &f);

void write_cross_dis_all(std::list<Organism> &organisms, std::ofstream &f);