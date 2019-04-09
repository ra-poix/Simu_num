#ifndef OPTION_H
#define OPTION_H


class Option{

    private:
        double strike;
        double horizon;

    public:
        double boundary_left();
        double boundary_right();
        double payoff(double x);

        double Horizon();
        double Strike();
        
};

#endif