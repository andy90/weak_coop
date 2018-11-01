#include "usefulfunctions.h"
#include "array_functions.h"
#include "constants.h"
#include "Organism.h"

double energy(std::vector<double> comp, std::vector<double> task) {
    double E = 0.0;
    double dissquare = normsquare(comp, task);
    if (dissquare<dis2square) {
        E = Emax - dissquare / D;
    }
    return E;
}

void calc_ntE_ntE2(CompCluster_Energy &CCE, std::vector<double> &ntE, std::vector<double> &ntE2) {

    // first clear the value of ntE and ntE2
    for (int i = 0; i < ntE.size(); i++) {
        ntE[i] = 0;
        ntE2[i] = 0;
    }

    for (auto i = CCE.begin(); i != CCE.end(); i++) {
        for (auto j = i->second.begin(); j != i->second.end(); j++) { // look at all the tasks the components are connected to
            ntE2[j->first] += 1;
            if (j->second > Estrong) {
                ntE[j->first] += 1;
            }
        }
    }


}

double costfunction(int N, std::vector<double> ntE, std::vector<double> ntE2) {
    double F = 0.0;
    if (ntE.size() != ntE2.size()) {
        std::cout << "size mismatch\n";
        return F;
    }

    int M = ntE.size();

    double t1 = 0.0; // this is the first term in the free energy functional
    for (int i = 0; i<M; i++) {
        if (ntE2[i]>0.0001) {
            t1 += 1;
        }
    }
    t1 = (M - t1); // shows how many edges are not connected


    double t2 = 0.0; // this is the second term in the free energy functional
    for (int i = 0; i<M; i++) {
        if (ntE[i]>0.0001) {
            t2 += 1;
        }
    }
    t2 = (M - t2); // shows how many edges are not strongly connected

    double t3 = N; // the third term of the functional

    F = 1.0*t1 + 1.0*t2 + 0.1*t3; // the total functional
    return F;

}

void add_comp(Comp_Task &components, CompCluster_Energy &CCE, std::vector<double> new_comp, int new_ic, std::map<int, double> energy_ic, coopflag cflag) {
    //Clusters clusters_connectedto;
    std::pair<Cluster, std::map<int, double>> all_cluster;
    //Cluster all; // stores all the cluster that is connected to new_ic
    all_cluster.second = energy_ic; // the total energy of the cluster "all", initialized with the energy of the new component
    all_cluster.first.push_back(new_ic); // put in the new added one first

    if (cflag == 1) { // only when the model is coop
        auto ii = CCE.begin();
        while (ii != CCE.end()) {
            //Cluster& cc = ii->first;
            int signal = 0;
            for (auto i = ii->first.begin(); i != ii->first.end(); i++) {
                int i_comp = *i;
                //std::vector<double> comp_i = components[i_comp];
                if (norm(new_comp, components[i_comp]) < dis_comp) {
                    for (auto j = ii->first.begin(); j != ii->first.end(); j++) { // push all the comp id in the cluster ii->first into all once it is connected to the new element
                        all_cluster.first.push_back(*j);
                    }
                    for (auto j = ii->second.begin(); j != ii->second.end(); j++) {
                        auto jtask = j->first;
                        auto jenergy = j->second;
                        if (all_cluster.second.find(jtask) != all_cluster.second.end()) {
                            all_cluster.second[jtask] += jenergy;
                        }
                        else {
                            all_cluster.second[jtask] = jenergy;
                        }
                    }
                    signal = 1;
                    break;
                }
            }

            if (signal != 0) {
                ii = CCE.erase(ii);
            }
            else {
                ii++;
            }
        }
    }

    CCE.push_back(all_cluster);

}


void delete_comp(CompCluster_Energy &CCE, Comp_info &tasks_for_comp, Comp_Task &components, int ic, coopflag cflag) {
    // ic is the component to be deleted
    // c_ic is the cluster ic belongs to, tasks for comp stores the task-comp energy, which is needed to reform
    // clusters
    Cluster c_ic;
    for (auto ii = CCE.begin(); ii != CCE.end(); ii++) {
        Cluster cc = ii->first;
        auto a = find(cc.begin(), cc.end(), ic);
        if (a != cc.end()) { // if we find ic in cc, record cc in c_ic. Then erase cc from CCE
            c_ic = cc;
            CCE.erase(ii);
            break;
        }
    }

    c_ic.remove(ic); // remove ic from the cluster c_ic

    CompCluster_Energy CCE_temp;
    for (auto i = c_ic.begin(); i != c_ic.end(); i++) {
        auto energy_i = tasks_for_comp[(*i)];
        add_comp(components, CCE_temp, components[(*i)], *i, energy_i, cflag); // re-add the not deleted components in the c_ic cluster 
    }

    for (auto i = CCE_temp.begin(); i != CCE_temp.end(); i++) {
        CCE.push_back(*i);
    }
    //CCE.merge(CCE_temp);
}


void generate_mutatelist(Rand_engine &generator, std::vector<int> &mlist) {
    // mlist is a vector of length of the number of components in a organism. now we try to see which components need to be mutates or duplicated
    int n = mlist.size();
    for (int i = 0; i < n; i++) {
        std::uniform_real_distribution<double> uni(0, 1);
        double dice = uni(generator); // use dice to determine which move to pick
        if (dice < mu / n) {
            mlist[i] = 1;
        }
        else if (dice < (mu + dup) / n) {
            mlist[i] = 2;
        }
        else if (dice < (mu + dup + loss) / n) {
            mlist[i] = 3;
        }
        else {
            mlist[i] = 0;
        }
    }
}

std::map<int, double> calc_tasksenergys(std::vector<double> &component, Comp_Task &tasks) {
    std::map<int, double> itasks_energys; // store the task the component is connected to
    int M = tasks.size();
    for (int j = 0; j<M; j++) {
        double energyij = energy(component, tasks[j]);
        if (energyij>0) {
            itasks_energys[j] = energyij; // component is is connected to task j, with energyij
        } // this is weak connection
    }

    return itasks_energys;
}

void mutate(std::list<Organism> &organisms, Rand_engine &generator) {
    std::list<Organism> orgas_new; // a temporary list to store the newly generated organisms
    for (auto i = organisms.begin(); i != organisms.end(); i++) {
        int ncomp = i->getN();
        int norg = i->getNorg(); // the number of organisms of species i
        for (int j = 0; j < norg; j++) {
            std::vector<int> mlist(ncomp);
            generate_mutatelist(generator, mlist);
            if (std::any_of(mlist.begin(), mlist.end(), [](int i) {return i != 0; })) {
                orgas_new.push_back(i->mutate_birth(generator, mlist));
            }
        }
    }

    for (auto i = orgas_new.begin(); i != orgas_new.end(); i++) {
        organisms.push_back(*i);
    }
}

void resample(std::list<Organism> &organisms, Rand_engine &generator) {
    // the following code makes selection and resample
    std::vector<double> fitness;
    for (auto i = organisms.begin(); i != organisms.end(); i++) {
        double cost = i->cost();
        double fit = exp(-beta * cost)*(i->getNorg()) / Norg;
        fitness.push_back(fit);
    }

    double fitness_av = 0;
    for (int i = 0; i < fitness.size(); i++) {
        fitness_av += fitness[i];
    }

    std::vector<double> p(fitness.size());
    for (int i = 0; i < fitness.size(); i++) {
        p[i] = fitness[i] / fitness_av;
    }

    std::discrete_distribution<int> dd(p.begin(), p.end()); // a discrete distribution based on the probability p
    std::vector<int> new_samplecount(fitness.size());
    for (int i = 0; i < Norg; i++) {
        int id = dd(generator);
        new_samplecount[id] += 1;
    }

    auto i_organism = organisms.begin();
    int i = 0;
    while (i_organism != organisms.end()) {
        i_organism->setNorg(new_samplecount[i]); // reset the number of organisms
        if (new_samplecount[i] == 0) {
            i_organism = organisms.erase(i_organism);
        }
        else {
            i_organism++;
        }
        i++;
    }
}

void mutateresample(std::list<Organism> &organisms, Rand_engine &generator, std::ofstream &f) {
    // the following code make the mutation
    std::list<Organism> orgas_new; // a temporary list to store the newly generated organisms
    for (auto i = organisms.begin(); i != organisms.end(); i++) {
        int ncomp = i->getN();
        int norg = i->getNorg(); // the number of organisms of species i
        for (int j = 0; j < norg; j++) {
            std::vector<int> mlist(ncomp);
            generate_mutatelist(generator, mlist);
            if (std::any_of(mlist.begin(), mlist.end(), [](int i) {return i != 0; })) {
                orgas_new.push_back(i->mutate_birth(generator, mlist));
            }
        }
    }

    for (auto i = orgas_new.begin(); i != orgas_new.end(); i++) {
        organisms.push_back(*i);
    }

    // the following code makes selection and resample
    std::vector<double> fitness;
    for (auto i = organisms.begin(); i != organisms.end(); i++) {
        double cost = i->cost();
        double fit = exp(-beta * cost)*(i->getNorg()) / Norg;
        fitness.push_back(fit);
    }

    double fitness_av = 0;
    for (int i = 0; i < fitness.size(); i++) {
        fitness_av += fitness[i];
    }

    std::vector<double> p(fitness.size());
    for (int i = 0; i < fitness.size(); i++) {
        p[i] = fitness[i] / fitness_av;
    }

    std::discrete_distribution<int> dd(p.begin(), p.end()); // a discrete distribution based on the probability p
    std::vector<int> new_samplecount(fitness.size());
    for (int i = 0; i < Norg; i++) {
        int id = dd(generator);
        new_samplecount[id] += 1;
    }

    auto i_organism = organisms.begin();
    int i = 0;
    while (i_organism != organisms.end()) {
        i_organism->setNorg(new_samplecount[i]); // reset the number of organisms
        if (new_samplecount[i] == 0) {
            i_organism = organisms.erase(i_organism);
        }
        else {
            i_organism++;
        }
        i++;
    }

    int n_success = orgas_new.size(); // the number of newly generated organisms
    int Nmutant_resampled = 0;
    int Nmutant_survived = 0;
    for (int i = (new_samplecount.size() - n_success); i < new_samplecount.size(); i++) {
        Nmutant_resampled += new_samplecount[i];
        if (new_samplecount[i] > 0) {
            Nmutant_survived += 1; // as long as the mutant survived, no matter how many children it has, count as 1.
        }
    }

    f << double(n_success) / Norg << " " << double(Nmutant_resampled) / Norg <<" "<<double(Nmutant_survived)/Norg<< "\n";
}

void write_completion_status(std::list<Organism> &organisms, std::ofstream &f) {
    
    std::vector<double> n_task(5); // stores the completion status of a task
    for (auto i_organism = organisms.begin(); i_organism != organisms.end(); i_organism++) {
        int M = i_organism->getM();
        std::vector<double> ntE(M), ntE2(M);
        calc_ntE_ntE2(i_organism->get_CCE(), ntE, ntE2);
        write_edge_percent(i_organism->get_CCE(), i_organism->get_tasks_for_comp(), ntE, ntE2, n_task, i_organism->getNorg());
    }

    for (int i = 0; i < 5; i++) {
        n_task[i] = n_task[i] / Norg;
        f << n_task[i] << " ";
    }
    f << "\n";
}

void write_edge_percent(CompCluster_Energy &CCE, Comp_info &tasks_for_comp, std::vector<double> &ntE, std::vector<double> &ntE2, std::vector<double> &n_task, int norg) {
    int n_nc = 0; // number of not completed tasks
    int n_spe_strong = 0; // strong specificly completed tasks
    int n_coop_strong = 0; // cooperatively specificly completed tasks
    int n_spe_coop = 0; // completed both cooply and strongly
    int n_weak = 0; // weakly completed tasks

    int M = ntE.size();
    std::vector<double> task_spe(M); // record the tasks that are strongly completed
    std::vector<double> task_coop(M); // record the tasks that are cooply completed


    for (auto i = CCE.begin(); i != CCE.end(); i++) {
        auto Ci = i->first;
        auto itasks_energys = i->second;
        for (auto j = itasks_energys.begin(); j != itasks_energys.end(); j++) {
            if (j->second > Estrong) { // see whether an element in Ci strongly completes a task. even a cluster is formed, if a component in the cluster strongly completes a task, it is counted as specificly completed
                if (find_if(Ci.begin(), Ci.end(), [&](int k)->bool {
                    if (tasks_for_comp[k][j->first] > Estrong) {
                        return 1;
                    }
                    else {
                        return 0;
                    }
                }) != Ci.end()) {
                    task_spe[j->first] += 1;
                }
                else {
                    task_coop[j->first] += 1;
                }
            }

        }
    }

    n_nc = M - count_if(ntE2.begin(), ntE2.end(), nonzero);
    n_weak = count_if(ntE2.begin(), ntE2.end(), nonzero) - count_if(ntE.begin(), ntE.end(), nonzero);
    n_spe_coop = count_if(task_spe.begin(), task_spe.end(), nonzero) + count_if(task_coop.begin(), task_coop.end(), nonzero) - count_if(ntE.begin(), ntE.end(), nonzero);
    n_spe_strong = count_if(task_spe.begin(), task_spe.end(), nonzero) - n_spe_coop;
    n_coop_strong = count_if(task_coop.begin(), task_coop.end(), nonzero) - n_spe_coop;


    n_task[0] += norg * n_nc; // stores the completion status of the task, weighted by the number of copies in an allel now. will be divided by Norg in the future
    n_task[1] += norg * n_weak;
    n_task[2] += norg * n_spe_coop;
    n_task[3] += norg * n_spe_strong;
    n_task[4] += norg * n_coop_strong;
}

bool nonzero(double a) {
    if (a > 0.1) {
        return 1;
    }
    else {
        return 0;
    }
}

void write_mutations(std::list<Organism> &organisms, std::ofstream &f) { // this is to write out the average mutation distance from the initial organism
    double n_mutations_av = 0;
    for (auto i_organism = organisms.begin(); i_organism != organisms.end(); i_organism++) {
        n_mutations_av += (i_organism->get_n_mutations())*(i_organism->getNorg()); // average mutations weighted by the the number of organisms
    }
    n_mutations_av = n_mutations_av / double(Norg);
    f << n_mutations_av << " \n";

}

std::vector<double> write_ntE_av(std::list<Organism> &organisms, int M, std::ofstream &f) {
    std::vector<double> ntE_av(M); // the average ntE at an instantaneous step
    for (auto i_organism = organisms.begin(); i_organism != organisms.end(); i_organism++) {
        std::vector<double> ntE(M), ntE2(M);
        calc_ntE_ntE2(i_organism->get_CCE(), ntE, ntE2);
        for (int i = 0; i < M; i++) {
            ntE_av[i] += ntE[i] * (i_organism->getNorg()) / double(Norg); // the ntE_av should be weighted correctly
        }

    }

    for (int i = 0; i < M; i++) {
        f << ntE_av[i] << " ";
    }
    f << "\n";
	return ntE_av;
}

void write_traj_mostprob(std::list<Organism> &organisms, std::ofstream &f, int isteps) {

    int norg_max = (organisms.begin())->getNorg(); // this is going to store the species with the maximum number of copies
    Organism the_organisim = *(organisms.begin());
    for (auto i_organism = organisms.begin(); i_organism != organisms.end(); i_organism++) {
        if (i_organism->getNorg() > norg_max) {
            norg_max = i_organism->getNorg();
            the_organisim = *i_organism;
        }
    }
    write_traj(the_organisim, f, isteps);
}

void write_traj(Organism &the_organism, std::ofstream &traj, int isteps) {
    int M = the_organism.getM();
    for (int i = 0; i<M; i++) {
        traj << isteps << " " << 0 << " ";
        Comp_Task tasks = the_organism.get_tasks();
        for (int j = 0; j<D; j++) {
            traj << tasks[i][j] << " ";
        }
        traj << "\n";
    }

    int l = 1;
    CompCluster_Energy CCE = the_organism.get_CCE();
    Comp_Task components = the_organism.get_components();
    for (auto i = CCE.begin(); i != CCE.end(); i++) {
        traj << isteps << " " << 1 << " ";
        for (auto j = i->first.begin(); j != i->first.end(); j++) {
            for (int k = 0; k < D; k++) {
                traj << components[*j][k] << " ";
            }
        }
        traj << "\n";
        l++;
    }
}

void writeN(std::list<Organism> &organisms, std::ofstream &f) {
    double Nav = 0; // the averaged number of components in the group of organism
    for (auto i_organism = organisms.begin(); i_organism != organisms.end(); i_organism++) {
        Nav += i_organism->getN()*(i_organism->getNorg());
    }
    Nav = Nav / double(Norg);
    f << Nav << " \n";
}

void write_cost_av(std::list<Organism> &organisms, std::ofstream &f2) {

    double cost_av = 0; // the average cost for this configuration
    double fitness_av = 0; // the average fitness for this configuration

    for (auto i_organism = organisms.begin(); i_organism != organisms.end(); i_organism++) {
        int M = i_organism->getM();
        std::vector<double> ntE(M), ntE2(M);
        int N = i_organism->getN();
        calc_ntE_ntE2(i_organism->get_CCE(), ntE, ntE2);
        cost_av += costfunction(N, ntE, ntE2)* (i_organism->getNorg()) / double(Norg);
        fitness_av += exp(-beta * costfunction(N, ntE, ntE2))* (i_organism->getNorg()) / double(Norg);
    }

    f2 << cost_av << " " << fitness_av << " ";
}

std::vector<double> generate_delta(double d, int D, Rand_engine &generator_fixed) {
    std::vector<double> delta_comp_ic(D); // the array to return, used to store the displacement between the tasks
    std::normal_distribution<double> norm_dis(0.0, 1.0);
    for (int j = 0; j<D; j++) {
        delta_comp_ic[j] = norm_dis(generator_fixed); // first generate the a vector of norm distribution
    }
    double dis = normsquare2(delta_comp_ic); // this is infact the squre of the distantce
    for (int j = 0; j<D; j++) {
        delta_comp_ic[j] = delta_comp_ic[j] / sqrt(dis)*d; // make sure the new distance is d
    }
    return delta_comp_ic;

}

void write_cluster_size(CompCluster_Energy &CCE, std::vector<double> &Ndis, int norg) {

	// a component might complete one task weakly, but it might collaborate with other components to complete another task
	// write out the cluster size forms by components
	int Nlimit = Ndis.size();//assume that the cluster size can not exceed Nlimit. Will put out the number of clusters with size from 0 to Nlimit-1
	for (auto i = CCE.begin(); i != CCE.end(); i++) {
		Cluster Ci = i->first;
		int isize = Ci.size();
		if (isize > Nlimit) {
			std::cout << "cluster size out of limit\n";
		}
		else {
			Ndis[isize - 1] += norg; // norg is the number of copies for this specific allel
		}
	}

}

void write_csize_av(std::list<Organism> &organisms, std::ofstream &f) {
	int Nlimit = 20;
	std::vector<double> Ndis(Nlimit);

	for (auto i_organism = organisms.begin(); i_organism != organisms.end(); i_organism++) {
		write_cluster_size(i_organism->get_CCE(), Ndis, i_organism->getNorg());
	}

	for (int i = 0; i < Nlimit; i++) {
		f << Ndis[i] / Norg << " ";
	}
	f << "\n";
}

void write_cross_dis(CompCluster_Energy &CCE, int norg, Comp_Task &tasks, std::ofstream &f) { // this writes out the cross reactivity between tasks and and distance between task and one component of the cluster
	for (auto i = CCE.begin(); i != CCE.end(); i++) {
		Cluster Ci = i->first;
		auto itask_energy = i->second;

		if ((Ci.size()>1) && (itask_energy.size() > 1)) {
			int n_strong = 0;
			std::vector<int> strongid;
			for (auto i = itask_energy.begin(); i != itask_energy.end(); i++) {
				if (i->second > Estrong) {
					n_strong += 1;
					strongid.push_back(i->first);
				}
			}



			double dis_task = 0; // the average distance between tasks that are being completed strongly
			double dis_tc = 0; // the average distance between tasks and components
			if (n_strong > 1) {
				// iterate over the strong tasks
				for (int i = 0; i < (n_strong - 1); i++) {
					for (int j = (i + 1); j < n_strong; j++) {
						dis_task += norm(tasks[strongid[i]], tasks[strongid[j]]);
					}
				}
				dis_task = dis_task / (n_strong*(n_strong - 1) / 2);
				f << tasks.size() << " " << Ci.size() << " " << norg << " " << n_strong << " " << dis_task << "\n";
			}



		}

	}


}

void write_cross_dis_all(std::list<Organism> &organisms, std::ofstream &f) {
	for (auto i_organism = organisms.begin(); i_organism != organisms.end(); i_organism++) {
		write_cross_dis(i_organism->get_CCE(), i_organism->getNorg(), i_organism->get_tasks(), f);
	}
}
