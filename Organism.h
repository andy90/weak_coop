#pragma once
#include "names.h"
#include <random>

class Organism
{
public:
    Organism(Comp_Task &components_init, Comp_Task &tasks_init, coopflag iflag_init, int norg_init); // the constructor of Organism. If we know the components of the Organism, the tasks its completing and its type (coop or not), we can get all the info we want 
    Organism(Comp_Task &components_init, CompCluster_Energy &CCE_init, Comp_info &tasks_for_comp_init, int norg_init, coopflag iflag_init, Comp_Task* p_tasks_init, int n_mutations_init);
    ~Organism();
    Organism mutate_birth(Rand_engine &generator, std::vector<int> &mlist); // give birth to organisms, where mutation, duplication and deletion of components might happen, after birth, the mother will die
    int getN(); // returns the number of components in this Organism
    void add_task(Comp_Task &tasks_new); // add a new task, all the private members of the class need to be changed
    int getNorg(); // returns the number of organisms of this species
    double cost(); // calculate the cost of the organism under the current fitness landscape
    void setNorg(int n); // reset the number of norg in this sepcies
    Comp_info& get_tasks_for_comp(); // return tasks for comp
    CompCluster_Energy& get_CCE(); // return CCE;
    int getM(); // get the number of tasks
    int get_n_mutations(); // get the number of mutations
    Comp_Task& get_tasks(); // get the tasks the organism is subject to
    Comp_Task& get_components(); // get the components
private:
    Comp_Task components;
    CompCluster_Energy CCE;
    Comp_info tasks_for_comp;
    int norg; // number of copies of this organism
    coopflag iflag;
    Comp_Task* p_tasks; // pointer to the task, should be updated when new tasks are introduced
    int n_mutations; // record the number of mutations from the initial organism
};

