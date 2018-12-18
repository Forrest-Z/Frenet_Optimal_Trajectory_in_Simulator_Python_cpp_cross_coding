// Transform from Frenet s,d coordinates to Cartesian x,y


/////////////////////////////////////////////////////////////////////////////
for (int i = -NUM_WAYPOINTS_BEHIND; i < NUM_WAYPOINTS_AHEAD; i++)
						coarse_waypoints_s.push_back(current_s); // index 0~9
						coarse_waypoints_x.push_back(map_waypoints_x[idx]);
						coarse_waypoints_y.push_back(map_waypoints_y[idx]);
						coarse_waypoints_dx.push_back(map_waypoints_dx[idx]);
						coarse_waypoints_dy.push_back(map_waypoints_dy[idx]);



////////////////////////////////////////////////////////////////////////////////////////////////////
interpolated_waypoints_x = interpolate_points(coarse_waypoints_s, coarse_waypoints_x, dist_inc, num_interpolation_points);//543					
					interpolated_waypoints_x = interpolate_points(coarse_waypoints_s, coarse_waypoints_x, dist_inc, num_interpolation_points);
					interpolated_waypoints_y = interpolate_points(coarse_waypoints_s, coarse_waypoints_y, dist_inc, num_interpolation_points);
					interpolated_waypoints_dx = interpolate_points(coarse_waypoints_s, coarse_waypoints_dx, dist_inc, num_interpolation_points);
					interpolated_waypoints_dy = interpolate_points(coarse_waypoints_s, coarse_waypoints_dy, dist_inc, num_interpolation_points);


///////////////////////////////////////////////////////////////////////////////////
					my_car.s    = pos_s;           // s position
					my_car.s_d  = s_dot;           // s dot - velocity in s
					my_car.s_dd = s_ddot;          // s dot-dot - acceleration in s
					my_car.d    = pos_d;           // d position
					my_car.d_d  = d_dot;           // d dot - velocity in d
					my_car.d_dd = d_ddot;          // d dot-dot - acceleration in d

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
for (string state: my_car.available_states)
vector<vector<double>> target_s_and_d = my_car.get_target_for_state(state, predictions, duration, car_just_ahead);
vector<vector<double>> possible_traj = my_car.generate_traj_for_target(target_s_and_d, duration);//JMT

//////////////////////////////////////////////////////////////////////
if (current_cost < best_cost) {
								best_cost = current_cost;
								best_frenet_traj = possible_traj;
								best_traj_state = state;
								best_target = target_s_and_d;
						}		
						
						
//////////////////////////////////////////////////////////////////////////////////////////////////////.
//coarse trajectory,
vector<double> coarse_s_traj, 
				coarse_x_traj, 
				coarse_y_traj,
				
				interpolated_s_traj,
				interpolated_x_traj, 
				interpolated_y_traj;





/////////////////////////////////////////////////////////////////////////////////////////////////////////
					msgJson["next_x"] = next_x_vals;
					msgJson["next_y"] = next_y_vals;










 
 ///////////////////////////////////////////////////////////////

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 