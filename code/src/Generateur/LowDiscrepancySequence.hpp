#ifndef LOWDISCREPANCYSEQUENCE_H
#define LOWDISCREPANCYSEQUENCE_H

#include <vector>
#include <algorithm>

class LowDiscrepancySequence{

    protected:
        long dim;

    public:
        LowDiscrepancySequence(long _dim): dim(_dim){};
        virtual std::vector<double> next(){
            return std::vector<double>();
        };

};

#endif