#include "Binary.hpp"

Binary::Binary(double _strike, double _horizon): 
    Option([_strike,_horizon] (double x) { 
                                if (x > _strike)
                                    return 1.0;
                                return 0.0;
            }),
    strike(_strike), horizon(_horizon){};
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