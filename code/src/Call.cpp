#include "Call.hpp"

using namespace std;


double Call::payoff(double x){
    return max(x-strike, 0.0);
}

double Call::boundary_left(double t, double x, double rate){
    return 0;
}

double Call::boundary_right(double t, double x, double rate){
    return x-strike*exp( -rate*(horizon-t) );
}