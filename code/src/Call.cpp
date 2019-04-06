#include "Call.hpp"

using namespace std;

double Call::payoff(double x){
    return max(x-strike, 0.0);
}