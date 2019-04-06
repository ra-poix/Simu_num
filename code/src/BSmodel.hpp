#include "model.hpp"

class BSmodel : public Model{

    public:
        BSmodel();
        BSmodel(double sigma, double r, double mu): sigma(sigma), r(r), mu(mu){}

        ~BSmodel();

        double Sigma();
        double Rate();

    private:
        double sigma;
        double r;

};