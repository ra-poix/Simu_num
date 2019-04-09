#ifndef MODEL_H
#define MODEL_H
class Model{

    private:

    public:
        Model();
        double Sigma(double t, double x);
        double Rate(double t, double x);
};

#endif