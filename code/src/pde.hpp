#ifndef PDE_H
#define PDE_H

#include <math.h>

#include "Option.hpp"
#include "Model.hpp"

class PDE{
    public:

        PDE(Option *_option, Model *_model): model(_model), option(_option){};
        ~PDE();

        double DV_coeff(double t, double x) const;
        double D2V_coeff(double t, double x) const;
        double V_coeff(double t, double x) const;
        double source_coeff(double t, double x) const;

        double boundary_left(double t, double x) const;
        double boundary_right(double t, double x) const;
        double init_cond(double x) const;

        Option* getOption();
        Model* getModel();

     private:
        Model *model;
        Option *option;
    
};

#endif