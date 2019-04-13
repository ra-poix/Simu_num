#include "option.hpp"

Option::Option(std::function<double (double) > _payoff): payoff(_payoff){}