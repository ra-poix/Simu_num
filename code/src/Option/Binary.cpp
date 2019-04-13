#include "Binary.hpp"

Binary::Binary(double _strike, double _horizon): 
    Option(_strike, _horizon, [_strike,_horizon] (double x) { 
                                if (x > _strike)
                                    return 1.0;
                                return 0.0;
            }) {};
double Binary::boundary_left(double t, double x, double rate){
    return 0;
}
double Binary::boundary_right(double t, double x, double rate){
    return 1;
}
