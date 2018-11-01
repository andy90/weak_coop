#include "Organism.h"
#include "usefulfunctions.h"
#include "array_functions.h"
#include "constants.h"
#include "names.h"
#include <vector>

Organism::Organism(Comp_Task &components_init, Comp_Task &tasks_init, coopflag iflag_init, int norg_init)
{
    components = components_init;
    iflag = iflag_init;
    norg = norg_init;
    p_tasks = &tasks_init;
    n_mutations = 0; // the initial organism is set to have 0 mutations
    // now construct CCE and tasks_for_comp
    for (auto i = components.begin(); i != components.end(); i++) {
        std::map<int, double> itasks_energys; // store the tasks and energys component i is connected to
        for (auto j = tasks_init.begin(); j != tasks_init.end(); j++) {
            double energyij = energy(i->second, j->second);
            if (energyij>0) {
                itasks_energys[j->first] = energyij; // component is is connected to task j, with energyij
            } // this is weak connection
        }
        tasks_for_comp[i->first] = itasks_energys; // itasks_energy can be empty.
    }

    // now construct CCE
    Comp_Task icomponents; // an intermediate strucutre for storing component coordinates
    for (auto i = components.begin(); i != components.end(); i++) {
        auto itasks_energys = tasks_for_comp[i->first]; // this is the most important information to get when we tries to add new components. We can define a function to get it
        add_comp(icomponents, CCE, components[i->first], i->first, itasks_energys, iflag); // icomponents is updated after the component is added. but I don't understand why I need to do this                                                                          // there is only 1 task initially, so does not matter whether the model is coop or not
        icomponents[i->first] = components[i->first];
    }
    
}

Organism::Organism(Comp_Task &components_init, CompCluster_Energy &CCE_init, Comp_info &tasks_for_comp_init, int norg_init, coopflag iflag_init, Comp_Task* p_tasks_init, int n_mutations_init) {
    components = components_init;
    CCE = CCE_init;
    tasks_for_comp = tasks_for_comp_init;
    norg = norg_init;
    iflag = iflag_init;
    p_tasks = p_tasks_init;
    n_mutations = n_mutations_init;
}

Organism::~Organism()
{
}

Organism Organism::mutate_birth(Rand_engine &generator, std::vector<int> &mlist)
{
    // mlist has too have non-zero elements, which means the children will have some mutations. it has the same number of elements as components
    CompCluster_Energy CCE_temp = CCE; // make a temporary copy of all the data and operate on them. notice this is a new copy which is independent of the original one
    Comp_info tasks_for_comp_temp = tasks_for_comp;
    Comp_Task components_temp = components;
    coopflag cflag = iflag;

    auto iter_comp = (components).begin(); // a iterator which points to the start of the components of this copy
    auto iter_taskcomp = (tasks_for_comp).begin(); // an iterator which points to the start of the task_for_complist of this copy
    for (int j = 0; j < components.size(); j++) { // try to mutate or duplicate the component in a copy one by one

        if (mlist[j] == 1) {
            // mutate the component
            int ic = iter_comp->first;
            std::lognormal_distribution<double> lg_dis(0.0, 1.0);
            double dis_move = lg_dis(generator)/len_scale; // the move distance for the mutation
            std::vector<double> delta_comp_ic1 = generate_delta(dis_move, D, generator); // the mutation in the coordinate of the component is decided by a normal distribution
            std::vector<double> comp_ic_new1(D); // stores the coordinate for the mutated component
            for (int i = 0; i < D; i++) {
                comp_ic_new1[i] = iter_comp->second[i] + delta_comp_ic1[i]; // the new coordinate for the component after the mutation
            }

            delete_comp(CCE_temp, tasks_for_comp_temp, components_temp, ic, cflag); // delete the component ic from CCE. don't need to update components  antd task_for_comp for now


            std::map<int, double> c1tasks_energys = calc_tasksenergys(comp_ic_new1, *p_tasks); // store the tasks and energys the new component ic is connected to			
            add_comp(components_temp, CCE_temp, comp_ic_new1, ic, c1tasks_energys, cflag); // add the mutated component ic. 
            components_temp[ic] = comp_ic_new1; // update the components list
            tasks_for_comp_temp[ic] = c1tasks_energys; //update tasks_for_comp list
        }
        else if (mlist[j] == 2) {
            // duplicate the component
            std::uniform_int_distribution<int> uni_int2(0, 200000000);
            int ic_toadd1 = uni_int2(generator); // this is the index of the child 1
            while (components_temp.find(ic_toadd1) != components_temp.end()) {
                ic_toadd1 = uni_int2(generator); // keep generating new one if it already exists in the components id list
            }

            add_comp(components_temp, CCE_temp, iter_comp->second, ic_toadd1, iter_taskcomp->second, cflag);
            components_temp[ic_toadd1] = iter_comp->second;
            tasks_for_comp_temp[ic_toadd1] = iter_taskcomp->second;
        }
        else if (mlist[j] == 3) {
            // delete the component
            int ic = iter_comp->first;
            delete_comp(CCE_temp, tasks_for_comp_temp, components_temp, ic, cflag); // delete the component ic from CCE.
            tasks_for_comp_temp.erase(ic); // delete the component from tasks_for_comp and components_temp
            components_temp.erase(ic);
        }

        advance(iter_comp, 1); // move the iterator to the next one
        advance(iter_taskcomp, 1);
    }
    norg -= 1; // delete one organism from this species

    Organism newborn(components_temp, CCE_temp, tasks_for_comp_temp, 1, iflag, p_tasks, n_mutations+1);
    return newborn;
}

int Organism::getN() {
    return components.size();
}

void Organism::add_task(Comp_Task &tasks) {
    int M = tasks.size();
    for (auto i = tasks_for_comp.begin(); i != tasks_for_comp.end(); i++) { // update tasks_for_comp. For the inner loop, just push in the newly added task is enough
        int id = i->first; // run through the component of each allel
        double energyij = energy(components[id], tasks[M - 1]); // tasks[M-1] is the newly added task. i_complist have the comp coordinates for the current allel
        if (energyij > 0) {
            i->second.insert(std::pair<int, double>(M - 1, energyij)); // component is is connected to task M-1, with energyij
        } // this is weak connection				
    }

    // update the CCE allel. the existing clusters can connect to the new task
    for (auto i = CCE.begin(); i != CCE.end(); i++) {
        Cluster Ci = i->first;
        for (auto j = Ci.begin(); j != Ci.end(); j++) { // run through each component contained by cluster Ci
            if (tasks_for_comp[*j].find(M - 1) != tasks_for_comp[*j].end()) {
                if (i->second.find(M - 1) != i->second.end()) {
                    i->second[M - 1] += tasks_for_comp[*j][M - 1];
                }
                else {
                    i->second[M - 1] = tasks_for_comp[*j][M - 1];
                }
            }
        }
    }

    p_tasks = &tasks; // update the pointer to the new tasks
}

int Organism::getNorg() {
    return norg;
}

double Organism::cost() {
    int M = p_tasks->size();
    std::vector<double> ntE(M), ntE2(M);
    calc_ntE_ntE2(CCE, ntE, ntE2); // calculate ntE and ntE2
    double cost = costfunction(components.size(), ntE, ntE2); // calculate cost
    return cost;
}

void Organism::setNorg(int n) {
    norg = n;
}

Comp_info& Organism::get_tasks_for_comp() {
    return tasks_for_comp;
}

CompCluster_Energy& Organism::get_CCE() {
    return CCE;
}

int Organism::getM() {
    return p_tasks->size();
}

int Organism::get_n_mutations() {
    return n_mutations;
}

Comp_Task& Organism::get_tasks() {
    return *p_tasks;
}

Comp_Task& Organism::get_components() {
    return components;
}