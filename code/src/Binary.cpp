#include "Binary.hpp"

double Binary::boundary_left(double t, double x, double rate){
    return 0;
}
double Binary::boundary_right(double t, double x, double rate){
    return 1;
}
double Binary::payoff(double x){
    if(x > strike)
        return 1;
    return 0;
}