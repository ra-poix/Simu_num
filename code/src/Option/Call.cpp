#include "Call.hpp"

#include <algorithm>

using namespace std;

Call::Call(double _strike, double _horizon): 
    Option([_strike,_horizon] (double x) { 
                                if (x > _strike)
                                    return x-_strike;
                                return 0.0;
            }),
    strike(_strike), horizon(_horizon){};

/*double Call::payoff(double x){
    return max(x-strike, 0.0);
}*/

double Call::boundary_left(double t, double x, double rate){
    return 0;
}

double Call::boundary_right(double t, double x, double rate){
    return x-strike*exp( -rate*(horizon-t) );
}