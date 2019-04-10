#include "fdm.hpp"

void FDM::calculate_step_sizes(){
    dx = x_max/static_cast<double>(J-1);
    dt = t_max/static_cast<double>(N-1);
}
void FDM::set_initial_conditions(){
                
    double cur_spot = 0.0;

    old_result.resize(J, 0.0);
    new_result.resize(J, 0.0);
    x_values.resize(J, 0.0);

    for (long pos=0; pos<J; pos++) {
        cur_spot = static_cast<double>(pos)*dx;
        old_result[pos] = pde->init_cond(cur_spot);
        x_values[pos] = cur_spot;
    }

    prev_t = 0.0;
    cur_t = 0.0;
}

void FDM::calculate_boundary_conditions(){
    new_result[0] = pde->boundary_left(prev_t, x_values[0]);
    new_result[J-1] = pde->boundary_right(prev_t, x_values[J-1]);
}

void FDM::step_march(){
    std::ofstream fdm_out("fdm.csv");

    while(cur_t < t_max) {
        cur_t = prev_t + dt;
        calculate_boundary_conditions();
        calculate_inner_domain();
        for (long pos=0; pos<J; pos++) {
          fdm_out << x_values[pos] << " " << prev_t << " " << new_result[pos] << std::endl;
        }
        
        old_result = new_result;
        prev_t = cur_t;
    }

    fdm_out.close();
}

FDMEuler::FDMEuler(double _x_max, double _t_max, long _J,  long _n, PDE *_pde): FDM(_x_max, _t_max, _J, _n, _pde) {
  calculate_step_sizes();
  set_initial_conditions();
}

void FDMEuler::calculate_inner_domain() {
  for (long pos=1; pos<J-1; pos++) {
    double dt_sig = dt * (pde->D2V_coeff(prev_t, x_values[pos]));
    double dt_sig_2 = dt * dx * 0.5 * (pde->DV_coeff(prev_t, x_values[pos]));

    alpha = dt_sig - dt_sig_2;
    beta = dx * dx - (2.0 * dt_sig) + (dt * dx * dx * (pde->V_coeff(prev_t, x_values[pos])));
    gamma = dt_sig + dt_sig_2;

    new_result[pos] = ( (alpha * old_result[pos-1]) + 
                      (beta * old_result[pos]) + 
                      (gamma * old_result[pos+1]) )/(dx*dx) - 
      (dt*(pde->source_coeff(prev_t, x_values[pos])));
  }
}

FDMCrankNik::FDMCrankNik(double _x_max, double _t_max, long _J,  long _n, PDE *_pde): FDM(_x_max, _t_max, _J, _n, _pde) {
  calculate_step_sizes();
  set_initial_conditions();
}

/**
 * 
 * 
 **/
void FDMCrankNik::calculate_inner_domain(){
  double sigma = Pde()->getModel()->Sigma(0.0,0.0);
  double mu = Pde()->getModel()->Rate(0.0,0.0) - sigma*sigma/2;
  double aux = (sigma*sigma*dt)/(2*dx*dx);
  //Calculer les coefs
  double qu = -aux/2 - (mu*dt)/(4*dx);
  double qm= 1 + aux;
  double qd = -aux/2 + (dt)/(4*dx);

  double pu = aux/2 + (mu*dt)/(4*dx);
  double pm = 1 - aux;
  double pd = aux/2 - (mu*dt)/(4*dx);
  //calculer terme de droite

  //Thomas algo pour terme de gauche
}
/*void FDMEuler::calculate_step_sizes(){
  dx = x_max/static_cast<double>(J-1);
  dt = t_max/static_cast<double>(N-1);
}*/

/*void FDMEuler::set_initial_conditions(){
  double cur_spot = 0.0;

  old_result.resize(J, 0.0);
  new_result.resize(J, 0.0);
  x_values.resize(J, 0.0);

  for (long pos=0; pos<J; pos++) {
    cur_spot = static_cast<double>(pos)*dx;
    old_result[pos] = pde->init_cond(cur_spot);
    x_values[pos] = cur_spot;
  }

  prev_t = 0.0;
  cur_t = 0.0;
}*/

/*void FDMEuler::calculate_boundary_conditions() {
  new_result[0] = pde->boundary_left(prev_t, x_values[0]);
  new_result[J-1] = pde->boundary_right(prev_t, x_values[J-1]);
}*/

/*void FDMEuler::step_march() { 
  std::ofstream fdm_out("fdm.csv");
  
  while(cur_t < t_max) {
    cur_t = prev_t + dt;
    calculate_boundary_conditions();
    calculate_inner_domain();
    for (long pos=0; pos<J; pos++) {
      fdm_out << x_values[pos] << " " << prev_t << " " << new_result[pos] << std::endl;
    }
    
    old_result = new_result;
    prev_t = cur_t;
  }

  fdm_out.close();
}*/