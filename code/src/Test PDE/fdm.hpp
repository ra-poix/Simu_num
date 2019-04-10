#ifndef FDM_H
#define FDM_H

#include "pde.hpp"
#include <vector>
#include <iostream>
#include <fstream>

class FDM {
    protected:
        PDE *pde;

        double x_max; 
        long J; 
        double dx; 
        std::vector<double> x_values; 

        double t_max; 
        long N; 
        double dt;

        double prev_t, cur_t;

        double alpha, beta, gamma;

        std::vector<double> new_result;
        std::vector<double> old_result;

        void calculate_step_sizes();
        void set_initial_conditions();
        void calculate_boundary_conditions();
        virtual void calculate_inner_domain() = 0;

    public:
        FDM(double _x_max, double _t_max, long _J,long _n, PDE *_pde)
            : x_max(_x_max), J(_J), t_max(_t_max), N(_n), pde(_pde) {};

        void step_march();
        PDE* Pde(){ return pde; };
};

class FDMEuler : public FDM{ //Euler explicite
    protected:
        void calculate_inner_domain();
    public:
        FDMEuler(double _x_max, double _t_max, long _J, long _n, PDE *_pde);
};

class FDMCrankNik : public FDM{ // Crank Nicholson
    protected:
        void calculate_inner_domain();
    public:
        FDMCrankNik(double _x_max, double _t_max, long _J, long _n, PDE *_pde);
};

#endif