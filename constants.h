#pragma once
#include <cmath>
// set the parameters, which are made global
const double sigma_task = 1.0; // define the length scale
const int D = 3;  // the dimension of the shape space
const double dis_theta2 = 0.05*sigma_task; // the criteria for weak activation
const double dis2square = dis_theta2 * dis_theta2; // the squre of the distheta2, not used after here for the 3D case
const double dis_theta = 0.2*dis_theta2;// the criteria for strong activation
const double Emax = dis_theta2 * dis_theta2 / D; // the maximum energy between the component and a task
const double Estrong = Emax - dis_theta * dis_theta / D; // the energy criterial for strong activation, there should be a coefficient in front of Emax, Estrong, and in the calculation of energy(.,.), but shouldn't really affect the final result
const double sigma_move = dis_theta; // the distance of monte carlo movement, or in other words, the mutation distance of the gene.
const double dis_comp = 1*dis_theta; // the distance betwen components to form a cluster
const double PI = 3.14159265358979323846;
const int Norg = 100; // this is the number of organisms
const double mu = 0.5;
const double dup = 0.05;
const double loss = 0.05;
const double beta = 1.0;
const double len_scale = exp(0.5) / sigma_move; // this parameter is used to scale the log_normal distribution
