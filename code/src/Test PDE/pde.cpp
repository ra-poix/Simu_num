#include "pde.hpp"

PDE::PDE(Option *_option, Model *_model): model(_model), option(_option){}

Option* PDE::getOption(){ return option; }
Model* PDE::getModel(){ return model; }

// Coefficient devant D2V
double PDE::D2V_coeff(double t, double x) const {
  double vol = model->Sigma(t,x);
  return 0.5*vol*vol*x*x; 
}

// Coefficient devant DV
double PDE::DV_coeff(double t, double x) const {
  return (model->Rate(t,x))*x;
}

// Coefficient devant V
double PDE::V_coeff(double t, double x) const {
  return -(model->Rate(t,x));
}

// Coefficient pour les dividendes
double PDE::source_coeff(double t, double x) const {
  return 0.0;
}

double PDE::boundary_left(double t, double x) const {
  return option->boundary_left(t,x, model->Rate(t,x));
}

double PDE::boundary_right(double t, double x) const {
  return option->boundary_right(t,x,model->Rate(t,x));
}

double PDE::init_cond(double x) const {
  return option->payoff(x);
}