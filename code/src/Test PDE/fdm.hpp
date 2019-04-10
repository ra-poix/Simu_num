#ifndef FDM_H
#define FDM_H

#include "pde.hpp"
#include <vector>
#include <iostream>

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

        virtual void calculate_step_sizes() = 0;
        virtual void set_initial_conditions() = 0;
        virtual void calculate_boundary_conditions() = 0;
        virtual void calculate_inner_domain() = 0;

    public:
        FDM(double _x_max, double _t_max, long _J,long _n, PDE *_pde)
            : x_max(_x_max), J(_J), t_max(_t_max), N(_n), pde(_pde) {};
        virtual void step_march() = 0;
};

class FDMEuler : public FDM{
    protected:
        void calculate_step_sizes();
        void set_initial_conditions();
        void calculate_boundary_conditions();
        void calculate_inner_domain();

    public:
        FDMEuler(double _x_max, double _t_max, long _J, long _n, PDE *_pde);
        void step_march();
};

#endif