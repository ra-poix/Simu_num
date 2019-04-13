#include "BSmodel.hpp"

BSmodel::BSmodel(double r,  double vol, double s): Model(  [r] (double t, double x) { return r;}, 
                                                            [vol] (double t, double x) { return vol;}, 
                                                            s) {
}