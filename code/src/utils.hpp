#include "Model/model.hpp"
#include "Option/option.hpp"
#include <cmath>


double mv_delta_bs(Option x, Model y, double *G , double e);
double mv_gamma_bs(Option x, Model y, double *G , double e);
double mv_vega_bs(Option x, Model y, double *G, double e );

double Tangent_delta_bs(Option x, Model y, double *G, double e );
double Tangent_vega_bs(Option x, Model y, double *G, double e );

double Tangent_delta_cev(Option x, Model y, double *G, double e );

double FD_delta_bs(Option o, Model m, double *x, double e);
double FD_vega_bs(Option o, Model m, double *x, double e);
double FD_gamma_bs(Option o, Model m, double *x, double e);

double FD_delta_cev(Option o, Model m, double *x, double e);
double FD_gamma_cev(Option o, Model m, double *x, double e);