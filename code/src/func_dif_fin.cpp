#include "Model/BSmodel.hpp"
#include "Option/Call.hpp"
#include "Option/Binary.hpp"
#include "MeanVar.hpp"

#include "Generateur/Generateur.hpp"
#include "Generateur/RandomGenerator.hpp"
#include "Generateur/LowDiscrepancyGenerator.hpp"
#include "Generateur/LowDiscrepancySequence.hpp"

#include <cmath>
#include <new>

//calculer payoff call a S
//Calculer payoff call a S+d
//return f(S+d) - f(S) / d

double mv_delta_bs(Option x, Model y, double G , double e){
 //   std::cout << y.rate(0.0,0.0) << "zeez";
    double S =  y.S()*exp((y.rate(0.0,0.0)-0.5*pow(y.sigma(0.0,0.0),2))*x.Horizon()+G*x.Horizon()*y.sigma(0.0,0.0));
    return exp(-y.rate(0.0,0.0)*x.Horizon()) * x.payoff(S) * G * x.Horizon() / (y.S() * y.sigma(0.0,0.0) * x.Horizon());
};

static double FD_delta_bs(Option o, Model m, double x, double e){
    double S0 = m.S();
    double sigma = m.sigma(0.0,0.0);
    double r = m.rate(0.0,0.0);


    double S1 = o.payoff( (S0+e) *  exp( (r - (1/2)*sigma*sigma) * o.Horizon() + sigma*x ) );
    double S2 = o.payoff( (S0-e) *  exp( (r - (1/2)*sigma*sigma) * o.Horizon() + sigma*x ) );
    return (S1 - S2)/(2*e);
}

static double FD_vega_bs(Option o, Model m, double x, double e){
    double S0 = m.S();
    double sigma = m.sigma(0.0,0.0) + e;
    double r = m.rate(0.0,0.0);

    double S1 = o.payoff( (S0) *  exp( (r - (1/2)*sigma*sigma) * o.Horizon() + (sigma)*x ) );
    sigma = m.sigma(0.0,0.0) - e;
    double S2 = o.payoff(  (S0) *  exp( (r - (1/2)*sigma*sigma) * o.Horizon() + sigma*x ) );
    return (S1 - S2)/(2*e);
}

static double FD_gamma_bs(Option o, Model m, double x, double e){
    double S0 = m.S();
    double sigma = m.sigma(0.0,0.0);
    double r = m.rate(0.0,0.0);
    printf("\nFD_gamma_bs s = %f vol = %f r = %f \n",S0,sigma, r);fflush(stdout);


    double S1 = o.payoff( (S0+e) *  exp( (r - (1/2)*sigma*sigma) * o.Horizon() + sigma*x ) );
    double S2 = o.payoff(  (S0-e) *  exp( (r - (1/2)*sigma*sigma) * o.Horizon() + sigma*x ) );
    double S3 = o.payoff( (S0) *  exp( (r - (1/2)*sigma*sigma) * o.Horizon() + sigma*x ) );
    std::cout << "\nFD_gamma_bs S1 / S2 / S3 :"<< S1 <<"/"<<S2 << "/" << S3 << std::endl;
    std::cout << "\nresult = " << (S1 + S2 - 2*S3) / (e*e) << std::endl;

    return (S1 + S2 - 2*S3) / (e*e);
}


//pointeur sur fonction
//Utilisation : function mon_pointeur[] = {maFonction1, maFonction 2, maFonction3, ...}
//type de retour : double
//params d'entrée : Option, Model, double, double
typedef double (*functions) (Option, Model, double, double);

bool enough_precise(MeanVar *mv, size_t size_mv) {
    for(int i = 0 ; i < size_mv; i++){
        if(!mv[i].is_enough_precise())
            return false;
    }
    return true;
}

void MonteCarlo(Option o, BSmodel m, Generateur g, functions *func, MeanVar **mv, size_t size_func){
    int compteur = 0;
    while(/*!enough_precise(*mv,size_func)*/compteur < 100000){
        double X = g.normale();
        for(int i = 0; i < size_func; i++){
            mv[i] -> maj( func[i](o,m,X,1) );        
        }
        compteur ++;
    }
};


int main(){

    functions f[]  = {mv_delta_bs};
    size_t size_f = 1;

    double S = 100;
    double r = 0.02;
    double vol = 0.3;
    BSmodel b(r, vol, S);

    double K = 100;
    double T = 1;
    Call c(K,T);

    double precision = 0.005;
    MeanVar *mv = new MeanVar[size_f];
    for(int i = 0; i < size_f; i++)
        mv[i] = MeanVar(precision);

    RandomGenerator g(5);

    MonteCarlo(c, b, g, f, &mv, size_f);
    std::cout << "\nFIN MC\n";

    std::cout << "mean = "<< mv[0].mean() << ", var = " << mv[0].var();
    return 0;
}


