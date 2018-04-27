#define MAX_STEER_ANGLE_IN_DEG 25
#define FIT_CURVE_ORDER 3
#define MPH_TO_MPS 0.44704
#define NUM_STATE_VAR 6 // [x y psi v cte epsi]
#define NUM_ACTUATOR_VAR 2

#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"

// for convenience
using json = nlohmann::json;
using namespace std;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.rfind("}]");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                        int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}

int main() {
  uWS::Hub h;

  // MPC is initialized here!
  MPC mpc;
  double latency = 0.1; //100ms : Latency to compensate actuator response time in ms. 
  const double Lf = 2.67;
  const double MAX_STEER_ANGE_IN_RADIAN = deg2rad(MAX_STEER_ANGLE_IN_DEG);

  h.onMessage([&mpc,&latency,&Lf,&MAX_STEER_ANGE_IN_RADIAN](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    string sdata = string(data).substr(0, length);
    // cout << sdata << endl;
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
		
        if (event == "telemetry") {
			
			
		  // cout<<"inside event routine"<<endl;
		  
		  
          // j[1] is the data JSON object
          vector<double> ptsx = j[1]["ptsx"];
          vector<double> ptsy = j[1]["ptsy"];
          double px = j[1]["x"];
          double py = j[1]["y"];
          double psi = j[1]["psi"];
          double v = j[1]["speed"];
		  double steer_prev = j[1]["steering_angle"];
		  double throttle_prev = j[1]["throttle"];

		  //Converting velocity to meters per second from mph
		  // v = v * MPH_TO_MPS;
		  
		  // cout<<"1"<<endl;
		  
          /*
          * TODO: Calculate steering angle and throttle using MPC.
          *
          * Both are in between [-1, 1].
          *
          */
		  
		  //1. Converting points to vehicles cordinate system
		  Eigen::VectorXd ptsx_cvt = Eigen::VectorXd::Zero(ptsx.size());
		  Eigen::VectorXd ptsy_cvt = Eigen::VectorXd::Zero(ptsy.size());
		  
		  for(size_t i=0;i<ptsx.size();i++){
			double cvtx = (ptsx[i] - px);
			double cvty = (ptsy[i] - py);
			ptsx_cvt[i] = cvtx*cos(-psi) - cvty*sin(-psi);
			ptsy_cvt[i] = cvtx*sin(-psi) + cvty*cos(-psi);
		  }
		  
		  
		  // cout<<"2"<<endl;
		  
		  //2. fit the curve using polyfit
		  auto coeffs = polyfit(ptsx_cvt,ptsy_cvt,FIT_CURVE_ORDER);
		  
		  // cout<<"3"<<endl;
		  
		  //3. calculate the cross track error and epsi
		  // double cte = polyeval(coeffs,px) - py;
		  //considering px , py , psi = 0 as we have converted the cordinate system to vehicle cordinate system
		  double cte = polyeval(coeffs,0);
		  
		  // cout<<"3.2"<<endl;
		  //Due to the sign starting at 0, the orientation error is -f'(x).
		  // double epsi = psi - atan(3*coeffs[3]*px*px + 2*coeffs[2]*px + coeffs[1]);
		  double epsi = - atan(coeffs[1]);
		  
		  Eigen::VectorXd state(NUM_STATE_VAR);
		  
		  // cout<<"4"<<endl;
		  
		  //4. Add Latency
		  // double new_x = px + v*cos(psi)*latency; // x = x0 + v*cos(psi)*latency
		  // double new_y = py + v*sin(osi)*latency; // y = y0 + v*sin(psi)*latency
		  // double new_psi = psi - (v/Lf)*steer_prev*latency; // psi = psi0 + (v/Lf)*delta*latency;
		  // double new_v = v + throttle_prev*latency; // v = v0 + a*latency;
		  
		  double new_x =  v*cos(psi)*latency; // x = x0 + v*cos(psi)*latency
		  double new_y =  v*sin(psi)*latency; // y = y0 + v*sin(psi)*latency
		  double new_psi = -(v/Lf)*steer_prev*latency; // psi = psi0 + (v/Lf)*delta*latency;
		  double new_v = v + throttle_prev*latency; // v = v0 + a*latency;
		  
		  // cout<<"before solve"<<endl;
		  
		  state << new_x,new_y,new_psi,new_v,cte,epsi;
		  // state << 0.0,0.0,0.0,v,cte,epsi;
		  
		  auto vars = mpc.Solve(state,coeffs);
		  
		  // cout<<"after solve"<<endl;
		  
          double steer_value = vars[0] / MAX_STEER_ANGE_IN_RADIAN;
          double throttle_value = vars[1];

          json msgJson;
          // NOTE: Remember to divide by deg2rad(25) before you send the steering value back.
          // Otherwise the values will be in between [-deg2rad(25), deg2rad(25] instead of [-1, 1].
          msgJson["steering_angle"] = steer_value;
          msgJson["throttle"] = throttle_value;

          //Display the MPC predicted trajectory 
          vector<double> mpc_x_vals;
          vector<double> mpc_y_vals;

          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Green line

		  for(int i=0;i<vars.size();i++){
			if(i%2==0){
				mpc_x_vals.push_back(vars[i]);
			}else{
				mpc_y_vals.push_back(vars[i]);
			}
		  }
		  
		  
          msgJson["mpc_x"] = mpc_x_vals;
          msgJson["mpc_y"] = mpc_y_vals;

          //Display the waypoints/reference line
          vector<double> next_x_vals;
          vector<double> next_y_vals;
		  
		  size_t num_points = 100;
		  int x_inc = 3;
		  
		  for(size_t i=0;i<num_points;i+=x_inc){
			next_x_vals.push_back(i);
			next_y_vals.push_back(polyeval(coeffs,i));
		  }
		  
          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Yellow line

          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;

          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          // std::cout << msg << std::endl;
          // Latency
          // The purpose is to mimic real driving conditions where
          // the car does actuate the commands instantly.
          //
          // Feel free to play around with this value but should be to drive
          // around the track with 100ms latency.
          //
          // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
          // SUBMITTING.
          this_thread::sleep_for(chrono::milliseconds(100));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
