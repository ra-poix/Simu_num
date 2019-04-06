#ifndef CALL_H
#define CALL_H

#include "option.hpp"
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>

class Call: public Option{

    public:
        double payoff(double x);

    private:
        double strike;

};

#endif