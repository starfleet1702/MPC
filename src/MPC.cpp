#define NUM_STATE_VAR 6 // [x y psi v cte epsi]
#define NUM_ACTUATOR_VAR 2

#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// DONE: Set the timestep length and duration
size_t N = 10;
double dt = 0.1;

/* ------------------START INDEX---------------------------
 * Variables to store starting index of different parameters
 * As it is required  to pass variables in vector form.
*/

size_t X_START = 0;
size_t Y_START = X_START + N;
size_t PSI_START = Y_START + N;
size_t V_START   = PSI_START + N;
size_t CTE_START = V_START + N;
size_t EPSI_START = CTE_START + N;
//FIXME
// size_t DELTA_START = EPSI_START + N - 1;
size_t DELTA_START = EPSI_START + N ;
size_t A_START = DELTA_START + N - 1;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

double ref_v = 8;//reference/desired velocity in mph(miles per hour)

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // DONE: implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
	
	// 1st value in fg is used to store cost value
	// COST = SE(CTE) + SE (EPSI) + SE(DIFF_VELOCITY) +
	//		  SQR(ABS_STEERING) + SQR(ABS_DELTA) +
	//        SQR(DIFF_IN_SEQ_STEERING) + SQR(DIFF_IN_SEQ_DELTA)
	//
	
	fg[0] = 0.0;
	
	//TODO : Add Weights
	
	// Adding SE of CTE , EPSI and DIFF_VELOCITY to cost function
	for(size_t t = 0; t < N; t++){
		fg[0] +=  0.41*CppAD::pow(vars[CTE_START+t],2);
		fg[0] +=  0.41*CppAD::pow(vars[EPSI_START+t],2);
		fg[0] +=  0.03*CppAD::pow(vars[V_START+t]-ref_v,2);
		
	}
	// Adding Square of value of delta and acceleration to minimize heavy use of actuators
	for(size_t t = 0; t < N-1; t++){
		fg[0] +=  7e-4*CppAD::pow(vars[DELTA_START+t],2);
		fg[0] +=  7e-4*CppAD::pow(vars[A_START+t],2);
		// try adding penalty for speed + steer
		fg[0] += 0.1*CppAD::pow(vars[DELTA_START + t] * vars[V_START+t], 2);
	}
	
	// Adding sequential difference of actuations for smooth driving
	for(size_t t = 1; t < N-1; t++){
		fg[0] +=  0.03*CppAD::pow(vars[DELTA_START+t] - vars[DELTA_START+t-1],2);
		fg[0] +=  0.0014*CppAD::pow(vars[A_START+t] - vars[A_START+t-1],2);
	}
	
	//---------------Vehicle Model Constrains----------------
	// It is Vehicle Model equations (Kinemetic 
	// model equations) which each point should follow 
	
	//Initial Constrains
	// We add 1 to each of the starting indices due to cost being located at
    // index 0 of `fg`. This bumps up the position of all the other values by 1.
	fg[1 + X_START] = vars[X_START];
	fg[1 + Y_START] = vars[Y_START];
	fg[1 + PSI_START] = vars[PSI_START];
	fg[1 + V_START] = vars[V_START];
	fg[1 + CTE_START] = vars[CTE_START];
	fg[1 + EPSI_START] = vars[EPSI_START];
	
	//Rest of the constraints
	for(size_t t = 1; t < N; t++){// starting at 1 , taking diff between 
	  //each step (current step) - prev step
		
	  // state at time t+1 (current)
      AD<double> x1 = vars[X_START + t];
	  AD<double> y1 = vars[Y_START + t];
      AD<double> psi1 = vars[PSI_START + t];
      AD<double> v1 = vars[V_START + t];
      AD<double> cte1 = vars[CTE_START + t];
      AD<double> epsi1 = vars[EPSI_START + t];
	  
	  // state at time t (prev)
      AD<double> x0 = vars[X_START + t - 1];
      AD<double> y0 = vars[Y_START + t - 1];
      AD<double> psi0 = vars[PSI_START + t - 1];
      AD<double> v0 = vars[V_START + t - 1];
      AD<double> cte0 = vars[CTE_START + t - 1];
      AD<double> epsi0 = vars[EPSI_START + t - 1];
	
	  // Only consider the actuation at time t (prev step).
	  AD<double> delta0 = vars[DELTA_START + t - 1];
	  AD<double> a0 = vars[A_START + t - 1];
	  
	  //degree 1 polyfit is a line , so only 2 coeff are there
	  AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * CppAD::pow(x0, 2)  + coeffs[3] * CppAD::pow(x0, 3);
	  
	  
	  //psides - will be slope of tanget to the curve -fit (path which we are following)
	  // m = tan(theta); theta = atan(m); m = df(x)/dx ; --> for line it will be coeffs[1]
	  // theta = atan(df(x)/dx);
	  AD<double> psides0 = CppAD::atan(coeffs[1] + 2 * coeffs[2] * x0 + 3 * coeffs[3] * CppAD::pow(x0, 2));
	  
	  //Rest of the model constraints
      fg[1 + X_START + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt); // xt = x0 + v*cos*dt
	  fg[1 + Y_START + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt); // yt = y0 + v*sin*dt
	  fg[1 + PSI_START + t] = psi1 - (psi0 - (v0/Lf)*delta0* dt);
	  fg[1 + V_START + t] = v1 - (v0 + a0* dt);
	  fg[1 + CTE_START + t] = cte1 - ((f0-y0) + (v0 * CppAD::sin(psi0) * dt));
	  fg[1 + EPSI_START + t] = epsi1 - ((psi0 - psides0) - v0 * delta0/Lf * dt);
	  
	}
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  
  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;
  
  //Extracting state variables
  double x = state[0];
  double y = state[1];
  double psi = state[2];
  double v = state[3];
  double cte = state[4];
  double epsi = state[5];

  // DONE : Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  size_t n_vars = NUM_STATE_VAR*N + NUM_ACTUATOR_VAR*(N-1);
  // DONE: Set the number of constraints
  size_t n_constraints = NUM_STATE_VAR * N;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }

  // Set the initial variable values
  vars[X_START] = x;
  vars[Y_START] = y;
  vars[PSI_START] = psi;
  vars[V_START] = v;
  vars[CTE_START] = cte;
  vars[EPSI_START] = epsi;
  
  
  //----------------LIMITS------------------
  
  // Lower and upper limits vector for variables
  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  
  // DONE : Set lower and upper limits for variables.
  
  // Setting all non-actuators upper and lowerlimits
  // to the max negative and positive values.
  for(size_t i=0;i<DELTA_START;i++){
	vars_lowerbound[i] = -1.0e19;
	vars_upperbound[i] = +1.0e19;
  }

  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).
  for (int i = DELTA_START; i < A_START; i++) {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332;
  }
  
  // Acceleration/decceleration upper and lower limits.
  for (int i = A_START; i < n_vars; i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }
  
  //----------------CONSTRAINS------------------
  
  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  
  //Initializing Constrains
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  //setting the constrain for initial state as this can't be changed
  constraints_lowerbound[X_START] = x;
  constraints_lowerbound[Y_START] = y;
  constraints_lowerbound[PSI_START] = psi;
  constraints_lowerbound[V_START] = v;
  constraints_lowerbound[CTE_START] = cte;
  constraints_lowerbound[EPSI_START] = epsi;

  constraints_upperbound[X_START] = x;
  constraints_upperbound[Y_START] = y;
  constraints_upperbound[PSI_START] = psi;
  constraints_upperbound[V_START] = v;
  constraints_upperbound[CTE_START] = cte;
  constraints_upperbound[EPSI_START] = epsi;
  
  
  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // cout<<"befor ipopt solve"<<endl;
  
  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= (solution.status == CppAD::ipopt::solve_result<Dvector>::success);

  // Cost
  auto cost = solution.obj_value;
  // std::cout << "Cost " << cost << std::endl;
  
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
  // return {solution.x[X_START + 1],   solution.x[Y_START + 1],
          // solution.x[PSI_START + 1], solution.x[V_START + 1],
          // solution.x[CTE_START + 1], solution.x[EPSI_START + 1],
		  // solution.x[DELTA_START],   solution.x[A_START]};
	
  // DONE: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  vector<double> output;
  
  output.push_back(solution.x[DELTA_START]);
  output.push_back(solution.x[A_START]);
  
  for(int i=1;i<N;i++){
	output.push_back(solution.x[X_START + i]);
	output.push_back(solution.x[Y_START	+ i]);
  }
  
  return output;
}
