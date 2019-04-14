#include "Model/BSmodel.hpp"
#include "Option/Call.hpp"
#include "Option/Binary.hpp"
#include "MeanVar.hpp"

#include "Generateur/Generateur.hpp"
#include "Generateur/RandomGenerator.hpp"
#include "Generateur/LowDiscrepancyGenerator.hpp"
#include "Generateur/LowDiscrepancySequence.hpp"

#include <cmath>
#include <iostream>
#include <fstream>
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
    return exp(-y.rate(0.0,0.0)*x.Horizon()) * x.payoff(S) * ( pow(G*sqrt(x.Horizon()),2) / (y.sigma(0.0,0.0)*x.Horizon()) -G*sqrt(x.Horizon()) - 1/y.sigma(0.0,0.0));
};


double Tangent_delta_bs(Option x, Model y, double G, double e ){
    double S = y.S()*exp((y.rate(0.0,0.0)-0.5*pow(y.sigma(0.0,0.0),2))*x.Horizon()+G*sqrt(x.Horizon())*y.sigma(0.0,0.0));
   // std::cout << S << std::endl;
    if(S>=x.Strike()) {return exp(-y.rate(0.0,0.0)*x.Horizon())S/y.S();} else {return 0;}
};

double Tangent_vega_bs(Option x, Model y, double G, double e ){
    double S = y.S()*exp((y.rate(0.0,0.0)-0.5*pow(y.sigma(0.0,0.0),2))*x.Horizon()+G*sqrt(x.Horizon())*y.sigma(0.0,0.0));
    if(S>=x.Strike()) {return exp(-y.rate(0.0,0.0)*x.Horizon())(G*sqrt(x.Horizon())-y.sigma(0.0,0.0)*x.Horizon())*S;} else {return 0;}
};



double Tangent_delta_cev(Option x, Model y, double G, double e ){
    euler scheme(x,y,gen,e);
    double S =G;
    if(S>=x.Strike()) {return exp(-y.rate(0.0,0.0)*x.Horizon())S/y.S();} else {return 0;}
};


static double FD_delta_bs(Option o, Model m, double x, double e){
    double S0 = m.S();
    double sigma = m.sigma(0.0,0.0);
    double r = m.rate(0.0,0.0);


    double S1 = o.payoff( (S0+e) *  exp( (r - 0.5*sigma*sigma) * o.Horizon() + sigma*x*sqrt(o.Horizon()) ) );
    double S2 = o.payoff( (S0-e) *  exp( (r - 0.5*sigma*sigma) * o.Horizon() + sigma*x*sqrt(o.Horizon()) ) );
    return (S1 - S2)/(2*e);
}

static double FD_vega_bs(Option o, Model m, double x, double e){
    double S0 = m.S();
    double sigma = m.sigma(0.0,0.0) + e;
    double r = m.rate(0.0,0.0);

    double S1 = o.payoff( (S0) *  exp( (r - 0.5*sigma*sigma) * o.Horizon() + (sigma)*x*sqrt(o.Horizon()) ) );
    sigma = m.sigma(0.0,0.0) - e;
    double S2 = o.payoff(  (S0) *  exp( (r - 0.5*sigma*sigma) * o.Horizon() + sigma*x*sqrt(o.Horizon()) ) );
    return (S1 - S2)/(2*e);
}

static double FD_gamma_bs(Option o, Model m, double x, double e){
    double S0 = m.S();
    double sigma = m.sigma(0.0,0.0);
    double r = m.rate(0.0,0.0);
    //printf("\nFD_gamma_bs x = %f S0 = %f vol = %f r = %f \n",x,S0,sigma, r);fflush(stdout);


    double S1 = o.payoff( (S0+e) *  exp( (r - 0.5*sigma*sigma) * o.Horizon() + sigma*x*sqrt(o.Horizon()) ) );
    double S2 = o.payoff(  (S0-e) *  exp( (r - 0.5*sigma*sigma) * o.Horizon() + sigma*x*sqrt(o.Horizon()) ) );
    double S3 = o.payoff( (S0) *  exp( (r - 0.5*sigma*sigma) * o.Horizon() + sigma*x*sqrt(o.Horizon()) ) );
   // std::cout << "\nFD_gamma_bs S1 / S2 / S3 :"<< S1 <<"/"<<S2 << "/" << S3 << std::endl;
   // std::cout << "\nresult = " << (S1 + S2 - 2*S3) / (e*e) << std::endl;

    return (S1 + S2 - 2*S3) / (e*e);
}


//pointeur sur fonction
//Utilisation : function mon_pointeur[] = {maFonction1, maFonction 2, maFonction3, ...}
//type de retour : double
//params d'entrée : Option, Model, double, double
typedef double (*functions) (Option, Model, double, double);
/*
bool enough_precise(MeanVar *mv) {
    for(int i = 0 ; i < size_mv; i++){
        if(!mv[i].is_enough_precise())
            return false;
    }
    
    return true;
}*/

void MonteCarlo(Option &o, BSmodel &m, RandomGenerator &g, functions func, MeanVar *mv){
    int compteur = 0;
    while(compteur < 1000000){
        compteur ++;
        if(compteur % 100000 == 0 && mv -> is_enough_precise())
            break;
        double X = g.normale();
        // printf("MonteCarlo g = %f \n",g);
        // std::cout << i << std::endl;fflush(stdout);
        mv -> maj( func(o,m,X,1) );
        mv -> maj( func(o,m,-X,1) );

    }

    std::cout << "compteur = " << compteur << std::endl;
};

void MonteCarlo(Option &o, CEVModel &m, RandomGenerator &g, functions func, MeanVar *mv){
    int compteur = 0;
    while(compteur < 1000000){
        compteur ++;
        if(compteur % 100000 == 0 && mv -> is_enough_precise())
            break;
        euler scheme(o,m,g);
        double X=scheme.solution(1000);
        // printf("MonteCarlo g = %f \n",g);
        // std::cout << i << std::endl;fflush(stdout);
        mv -> maj( func(o,m,X,1) );
        //mv -> maj( func(o,m,-X,1) );
        
    }
    
    std::cout << "compteur = " << compteur << std::endl;
};

int main(){

    functions f  = Tangent_vega_bs;
    int size_f = 1;

    double S;
    double r = 0.02;
    double vol = 0.4;
    BSmodel b(r, vol, S);

    double K = 100;
    double T = 1;
    Call c(K,T);

    double precision = 0.001;
    MeanVar *mv = new MeanVar;
    *mv = MeanVar(precision);

    RandomGenerator g(time(NULL));
    std::ofstream fdm_out("fdm.csv");
    fdm_out << "S" << ";" << "Y" << std::endl;
    for(S = 10 ; S < 200 ; S++){
        BSmodel b(r,vol,S);
        MonteCarlo(c, b, g, f, mv);
        std::cout<< "S = " << S << " mean = "<< mv->mean() << ", var = " << mv -> var() << std::endl;
        fdm_out << S << ";" << mv->mean() << std::endl;
        *mv = MeanVar(precision);
    }

    delete mv;

    //RandomGenerator g(time(NULL));

    std::cout << "\nFIN MC\n";
    std::cout << "mean = "<< mv[0].mean() << ", var = " << mv[0].var() << std::endl;
  //  std::cout << "mean = "<< mv[1].mean() << ", var = " << mv[1].var();


    return 0;
}


