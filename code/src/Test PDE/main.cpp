#include "fdm.hpp"
#include "pde.hpp"
#include "Call.hpp"
#include "BSmodel.hpp"

int main(){

    double s = 100;
    double T=1;
    Call *c =  new Call(s,T); // call strike 100 horizon 1an

    double vol = 0.3;
    double r = 0.02;
    double s0 = 100;
    BSmodel *mod = new BSmodel(vol,r,s0);

    PDE *pde = new PDE(c,mod);

    long x_max = 1000;
    long div_x = 1000;
    long div_t = 1000;
    FDMEuler fdm(x_max, c->Horizon(), div_x, div_t, pde);

    fdm.step_march();
    
    return 0;
}