#ifndef LOWDISCREPANCYSEQUENCE_H
#define LOWDISCREPANCYSEQUENCE_H

#include <vector>
#include <algorithm>
#include <functional>

class LowDiscrepancySequence{

    protected:
        long dim;

    public:
        LowDiscrepancySequence(std::function<std::vector<double> ()> _next ,long _dim): next(_next) ,dim(_dim){};

        std::function<std::vector<double> ()> next;

};

#endif