/**
 * Contient toutes les fonctions utilisées pour simuler les grecques
 * Malliavin, Différence finie, Processus tangent 
 */ 

#include "utils.hpp"


double mv_delta_bs(Option x, Model y, double *G , double e){
    double rate = y.rate(0.0,0.0);
    double sigma = y.sigma(0.0,0.0);
    double T = x.Horizon();
    double g = *G;
    double S = y.S()*exp((rate-0.5*pow(sigma,2))*T+g*sqrt(T)*sigma);
    return exp(-rate*T)*x.payoff(S)*g*sqrt(T)/(y.S()*sigma*x.Horizon());
};

double mv_gamma_bs(Option x, Model y, double *G , double e){
    double rate = y.rate(0.0,0.0);
    double sigma = y.sigma(0.0,0.0);
    double T = x.Horizon();
    double g = *G;

    double S = y.S()*exp((rate-0.5*pow(sigma,2))*T+g*sqrt(T)*sigma);

    double d1 = 0.5*(x.payoff(S+e)-x.payoff(S-e))/e;
    double d2 = 0.5*(x.payoff(y.S()*exp(rate*T) +e)-x.payoff(y.S()*exp(rate*T) -e))/e ;

    double num = (d1 * S - x.payoff(S) - y.S() * exp(rate*T)* d2 + x.payoff(y.S()*exp(rate*T)) ) *g;
    double den = (y.S() * y.S() * sigma * T );

    //return exp(-rate*T)*x.payoff(S)*(pow(g*sqrt(T),2)/(sigma*T)-g*sqrt(T)-1/sigma)/(pow(y.S(),2)*sigma*T);
    
    return num/den;
};

double mv_gamma_bs2(Option x, Model y, double *G , double e){
    double rate = y.rate(0.0,0.0);
    double sigma = y.sigma(0.0,0.0);
    double T = x.Horizon();
    double g = *G;

    double S = y.S()*exp((rate-0.5*pow(sigma,2))*T+g*sqrt(T)*sigma);

    double d1 = 0.5*(x.payoff(S+e)-x.payoff(S-e))/e;
    double d2 = 0.5*(x.payoff(y.S()*exp(rate*T) +e)-x.payoff(y.S()*exp(rate*T) -e))/e ;

    double num = (d1 * S - x.payoff(S) - y.S() * exp(rate*T)* d2 + x.payoff(y.S()*exp(rate*T)) * g);
    double den = (y.S() * y.S() * sigma * T );

    return exp(-rate*T)*x.payoff(S)*(pow(g*sqrt(T),2)/(sigma*T)-g*sqrt(T)-1/sigma)/(pow(y.S(),2)*sigma*T);
    
    //return num/den;
};

double mv_vega_bs(Option x, Model y, double *G, double e ){
    double g = *G;

    double S = y.S()*exp((y.rate(0.0,0.0)-0.5*pow(y.sigma(0.0,0.0),2))*x.Horizon()+g*sqrt(x.Horizon())*y.sigma(0.0,0.0));
    return exp(-y.rate(0.0,0.0)*x.Horizon()) * x.payoff(S) * ( pow(g*sqrt(x.Horizon()),2) / (y.sigma(0.0,0.0)*x.Horizon()) -g*sqrt(x.Horizon()) - 1/y.sigma(0.0,0.0));
};


double Tangent_delta_bs(Option x, Model y, double *G, double e ){
    double g = *G;
    double S = y.S()*exp((y.rate(0.0,0.0)-0.5*pow(y.sigma(0.0,0.0),2))*x.Horizon()+g*sqrt(x.Horizon())*y.sigma(0.0,0.0));
    if(S>=x.Strike()) {
        return exp(-y.rate(0.0,0.0)*x.Horizon())*S/y.S();
    } else {
        return 0;
    }
};

double Tangent_vega_bs(Option x, Model y, double *G, double e ){
    double g = *G;
    double S = y.S()*exp((y.rate(0.0,0.0)-0.5*pow(y.sigma(0.0,0.0),2))*x.Horizon()+g*sqrt(x.Horizon())*y.sigma(0.0,0.0));
    if(S>=x.Strike()) {
        return exp(-y.rate(0.0,0.0)*x.Horizon()) * (g*sqrt(x.Horizon())-y.sigma(0.0,0.0)*x.Horizon())*S;
    }else 
       return 0;
};

double Tangent_delta_cev(Option x, Model y, double *G, double e ){
    double S = *G;
    if(S>=x.Strike()) {
        return exp(-y.rate(0.0,0.0)*x.Horizon())*S/y.S();
    } else {
        return 0;
    }
};


double FD_delta_bs(Option o, Model m, double *x, double e){
    double S0 = m.S();
    double sigma = m.sigma(0.0,0.0);
    double r = m.rate(0.0,0.0);

    double S1 = o.payoff( (S0+e) *  exp( (r - 0.5*sigma*sigma) * o.Horizon() + sigma*(*x)*sqrt(o.Horizon()) ) );
    double S2 = o.payoff( (S0-e) *  exp( (r - 0.5*sigma*sigma) * o.Horizon() + sigma*(*x)*sqrt(o.Horizon()) ) );
    return (S1 - S2)/(2*e)* exp(-r*o.Horizon());
}

double FD_vega_bs(Option o, Model m, double *x, double e){
    double S0 = m.S();
    double sigma1 = m.sigma(0.0,0.0) + e;
    double sigma2=m.sigma(0.0,0.0) -e;
    double r = m.rate(0.0,0.0);

    double S1 = o.payoff( (S0) *  exp( (r - 0.5*sigma1*sigma1) * o.Horizon() + (sigma1)*(*x)*sqrt(o.Horizon()) ) );
    double S2 = o.payoff(  (S0) *  exp( (r - 0.5*sigma2*sigma2) * o.Horizon() + sigma2*(*x)*sqrt(o.Horizon()) ) );
    return (S1 - S2)/(2*e)* exp(-r*o.Horizon());
}

double FD_gamma_bs(Option o, Model m, double *x, double e){
    double S0 = m.S();
    double sigma = m.sigma(0.0,0.0);
    double r = m.rate(0.0,0.0);

    double S1 = o.payoff( (S0+e) *  exp( (r - 0.5*sigma*sigma) * o.Horizon() + sigma*(*x)*sqrt(o.Horizon()) ) );
    double S2 = o.payoff(  (S0-e) *  exp( (r - 0.5*sigma*sigma) * o.Horizon() + sigma*(*x)*sqrt(o.Horizon()) ) );
    double S3 = o.payoff( (S0) *  exp( (r - 0.5*sigma*sigma) * o.Horizon() + sigma*(*x)*sqrt(o.Horizon()) ) );

    return (S1 + S2 - 2*S3) / (e*e) * exp(-r*o.Horizon());
}

double FD_gamma_cev(Option o, Model m, double *x, double e){
    double S0 = m.S();
    double S1 = o.payoff( *(x+1) );
    double S2 = o.payoff( *(x+2) );
    double S3 = o.payoff(*x);

    return (S1 + S2 - 2*S3) / (e*e)* exp(-m.rate(0.0,0.0)*o.Horizon());
}

//ERREUR DANS LA FORMULE
double FD_delta_cev(Option o, Model m, double *x, double e){
    double S0 = m.S();

    double S1 = o.payoff( *(x+1) );
    double S2 = o.payoff( *(x+2));
    return (S1 - S2)/(2*e)* exp(-m.rate(0.0,0.0)*o.Horizon());
} 
