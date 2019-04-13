#include "option.hpp"

Option::Option(double _horizon, double _strike,std::function<double (double) > _payoff): 
                                            payoff(_payoff), 
                                            strike(_strike), 
                                            horizon(_horizon){}