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
  
  std::vector<double> right_term = calculate_right_term();
  
  std::vector<double> c, b, a;
  
  c.resize(J,0.0);
  b.resize(J,0.0);
  a.resize(J-1,0.0);
  double sigma, mu, aux;
  for(int i = 0; i < J; i++){
    sigma = Pde()->getModel()->sigma(cur_t,x_values[i]);
    mu = Pde()->getModel()->rate(cur_t,x_values[i]) - sigma*sigma/2;
    aux = (sigma*sigma*dt)/(2*dx*dx);

    b[i] = 1 + aux;

    if(i==J-1){
      a[J-2] = -aux/2 + (mu*dt)/(4*dx);
    }else if(i==0){
      c[i] = -aux/2 - (mu*dt)/(4*dx);
    }else{
      a[i-1] = -aux/2 + (mu*dt)/(4*dx);
      c[i] = -aux/2 - (mu*dt)/(4*dx);
    }
  }

  c[0] = c[0] / b[0];
  old_result[0] = old_result[0] / b[0];
  for(int i=1; i < J; i++){
    c[i] = c[i] / (b[i] - a[i-1]*c[i-1]);
    old_result[i] = (old_result[i] - a[i-1] * old_result[i-1])  / (b[i] - a[i-1]*c[i-1]);
  }
  new_result[J-1] = old_result[J-1];
  for(int i=J-2; i >= 0; i--){
    new_result[i] = old_result[i] - c[i]*new_result[i+1];
  }
}

std::vector<double> FDMCrankNik::calculate_right_term(){
  
  std::vector<double> right_term;
  right_term.resize(J,0.0);
  right_term[0] = 1;
  right_term[J-1] = 1;

  double sigma, mu, aux, pu, pm, pd;
  sigma = Pde()->getModel()->sigma(cur_t,0);
  mu = Pde()->getModel()->rate(cur_t,0) - sigma*sigma/2;
  aux = (sigma*sigma*dt)/(2*dx*dx);

  pu = aux/2 + (mu*dt)/(4*dx);
  pm = 1 - aux;
  pd = aux/2 - (mu*dt)/(4*dx);
  right_term[0] = pd*old_result[0] + pm*old_result[0] + pu*old_result[1];

  sigma = Pde()->getModel()->sigma(cur_t,x_values[J-1]);
  mu = Pde()->getModel()->rate(cur_t,x_values[J-1]) - sigma*sigma/2;
  aux = (sigma*sigma*dt)/(2*dx*dx);

  pu = aux/2 + (mu*dt)/(4*dx);
  pm = 1 - aux;
  pd = aux/2 - (mu*dt)/(4*dx);
  right_term[J-1] = pd*old_result[J-2] + pm*old_result[J-1] + pu*old_result[J-1];
  
  for(int i = 1; i < J-1; i++){
    sigma = Pde()->getModel()->sigma(cur_t,i*dx);
    mu = Pde()->getModel()->rate(cur_t,i*dx) - sigma*sigma/2;
    aux = (sigma*sigma*dt)/(2*dx*dx);

    pu = aux/2 + (mu*dt)/(4*dx);
    pm = 1 - aux;
    pd = aux/2 - (mu*dt)/(4*dx);

    right_term[i] = pd*old_result[i-1] + pm*old_result[i] + pu*old_result[i+1];
  }
  return right_term;
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