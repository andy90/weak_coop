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


typedef std::map<int, std::vector<double> > Comp_Task; // stores the id and coordinates for tasks and components
typedef std::list< std::pair< std::list<int>, std::map<int, double> > > CompCluster_Energy; // stores the clusters of components, which tasks each cluster is connected to and the energy of connection
typedef std::list<int> Cluster; // cluster formed by components
typedef std::list<std::list<int> > Clusters;
typedef std::map<int, std::map<int, double> > Comp_info; // information about each component connects to which task with which energy

typedef std::list<  Comp_Task  > Comp_Task_list; // stores the of coordinate of the components for the list of organisms, int is the id of the organism 
typedef std::list<  CompCluster_Energy  > CCE_list; // stores the CCE for the list of organisms, int is the id of the organism
typedef std::list<  Comp_info > Comp_info_list; // stores the how components connect to the task for the list of organisms.
typedef std::list<  int > number_list; // stores the number of each allel for the list of organisms
typedef std::list<int> coopflag_list; // a 0 or 1 flag which identifies whether a organism is coop or not
typedef int coopflag;

typedef std::mt19937_64 Rand_engine; // define a name for the random engine