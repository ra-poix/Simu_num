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

double mv_delta_bs(Option x, Model y, double G , double e){
    double S = y.S()*exp((y.rate(0.0,0.0)-0.5*pow(y.sigma(0.0,0.0),2))*x.Horizon()+G*sqrt(x.Horizon())*y.sigma(0.0,0.0));

   // std::cout << " g= "<< G << " horizon = " << S   << std::endl;
    return exp(-y.rate(0.0,0.0)*x.Horizon())*x.payoff(S)*G*sqrt(x.Horizon())/(y.S()*y.sigma(0.0,0.0)*x.Horizon());
};

double mv_gamma_bs(Option x, Model y, double G , double e){
    double S = y.S()*exp((y.rate(0.0,0.0)-0.5*pow(y.sigma(0.0,0.0),2))*x.Horizon()+G*sqrt(x.Horizon())*y.sigma(0.0,0.0));
    return exp(-y.rate(0.0,0.0)*x.Horizon())*x.payoff(S)*(pow(G*sqrt(x.Horizon()),2)/(y.sigma(0.0,0.0)*x.Horizon())-G*sqrt(x.Horizon())-1/y.sigma(0.0,0.0))/(pow(y.S(),2)*y.sigma(0.0,0.0)*x.Horizon());
};

double mv_vega_bs(Option x, Model y, double G, double e ){
    double S = y.S()*exp((y.rate(0.0,0.0)-0.5*pow(y.sigma(0.0,0.0),2))*x.Horizon()+G*sqrt(x.Horizon())*y.sigma(0.0,0.0));
    return exp(-y.rate(0.0,0.0)*x.Horizon())*x.payoff(S)*(pow(G*sqrt(x.Horizon()),2)/(y.sigma(0.0,0.0)*x.Horizon())-G*sqrt(x.Horizon())-1/y.sigma(0.0,0.0));
};


double Tangent_delta_bs(Option x, Model y, double G, double e ){
    double S = y.S()*exp((y.rate(0.0,0.0)-0.5*pow(y.sigma(0.0,0.0),2))*x.Horizon()+G*sqrt(x.Horizon())*y.sigma(0.0,0.0));
    if(S>=x.Strike()) {return S/y.S();} else {return 0;}
};

double Tangent_vega_bs(Option x, Model y, double G, double e ){
    double S = y.S()*exp((y.rate(0.0,0.0)-0.5*pow(y.sigma(0.0,0.0),2))*x.Horizon()+G*sqrt(x.Horizon())*y.sigma(0.0,0.0));
    if(S>=x.Strike()) {return (G*sqrt(x.Horizon())-y.sigma(0.0,0.0)*x.Horizon())*S/y.S();} else {return 0;}
};

/*double Tangent_delta_cev(Option x, Model y, double G, double e ){
    euler scheme(x,y,gen,e);
    double S = y.S()*exp((y.rate(0.0,0.0)-0.5*pow(y.sigma(0.0,0.0),2))*x.Horizon()+G*sqrt(x.Horizon())*y.sigma(0.0,0.0));
    if(S>=x.Strike()) {return S/y.S();} else {return 0;}
};*/ //A TRAVAILLER !!!!!


static double FD_delta_bs(Option o, Model m, double x, double e){
    double S0 = m.S();
    double sigma = m.sigma(0.0,0.0);
    double r = m.rate(0.0,0.0);


    double S1 = o.payoff( (S0+e) *  exp( (r - 0.5*sigma*sigma) * o.Horizon() + sigma*x ) );
    double S2 = o.payoff( (S0-e) *  exp( (r - 0.5*sigma*sigma) * o.Horizon() + sigma*x ) );
    return (S1 - S2)/(2*e);
}

static double FD_vega_bs(Option o, Model m, double x, double e){
    double S0 = m.S();
    double sigma = m.sigma(0.0,0.0) + e;
    double r = m.rate(0.0,0.0);

    double S1 = o.payoff( (S0) *  exp( (r - 0.5*sigma*sigma) * o.Horizon() + (sigma)*x ) );
    sigma = m.sigma(0.0,0.0) - e;
    double S2 = o.payoff(  (S0) *  exp( (r - 0.5*sigma*sigma) * o.Horizon() + sigma*x ) );
    return (S1 - S2)/(2*e);
}

static double FD_gamma_bs(Option o, Model m, double x, double e){
    double S0 = m.S();
    double sigma = m.sigma(0.0,0.0);
    double r = m.rate(0.0,0.0);
    //printf("\nFD_gamma_bs x = %f S0 = %f vol = %f r = %f \n",x,S0,sigma, r);fflush(stdout);


    double S1 = o.payoff( (S0+e) *  exp( (r - 0.5*sigma*sigma) * o.Horizon() + sigma*x ) );
    double S2 = o.payoff(  (S0-e) *  exp( (r - 0.5*sigma*sigma) * o.Horizon() + sigma*x ) );
    double S3 = o.payoff( (S0) *  exp( (r - 0.5*sigma*sigma) * o.Horizon() + sigma*x ) );
   // std::cout << "\nFD_gamma_bs S1 / S2 / S3 :"<< S1 <<"/"<<S2 << "/" << S3 << std::endl;
   // std::cout << "\nresult = " << (S1 + S2 - 2*S3) / (e*e) << std::endl;

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

void MonteCarlo(Option o, BSmodel m, RandomGenerator g, functions *func, MeanVar **mv, int size_func){
    int compteur = 0;
    while(compteur < 100000){
        compteur ++;
        if(compteur % 1000 == 0 && enough_precise(*mv,size_func))
            break;
        double X = g.normale();
       // printf("MonteCarlo g = %f \n",g);
        for(int i = 0; i < size_func; i++){
           // std::cout << i << std::endl;fflush(stdout);
            mv[i] -> maj( func[i](o,m,X,1) );
            mv[i] -> maj( func[i](o,m,-X,1) );
        }
        
    }

    std::cout << "compteur = " << compteur << std::endl;
};


int main(){

    functions f[]  = {FD_delta_bs};
    int size_f = 1;

    double S;
    double r = 0.02;
    double vol = 0.4;
    BSmodel b(r, vol, S);

    double K = 100;
    double T = 10;
    Call c(K,T);

    double precision = 0.001;
    MeanVar *mv = new MeanVar[size_f];
    for(int i = 0; i < size_f; i++)
        mv[i] = MeanVar(precision);
    RandomGenerator g(time(NULL));
    for(S = 50 ; S < 150 ; S++){
        BSmodel b(r,vol,S);
        MonteCarlo(c, b, g, f, &mv, size_f);
        std::cout << "mean = "<< mv[0].mean() << ", var = " << mv[0].var() << std::endl;
    }

    //RandomGenerator g(time(NULL));

    

    std::cout << "\nFIN MC\n";
    std::cout << "mean = "<< mv[0].mean() << ", var = " << mv[0].var() << std::endl;
  //  std::cout << "mean = "<< mv[1].mean() << ", var = " << mv[1].var();


    return 0;
}


