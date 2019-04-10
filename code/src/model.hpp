#ifndef MODEL_H
#define MODEL_H
class Model{

    public:
        virtual double Sigma(double t, double x) = 0;
        virtual double Rate(double t, double x) = 0;
};

#endif