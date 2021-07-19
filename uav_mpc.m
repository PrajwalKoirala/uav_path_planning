clear all; close all;

addpath('/MATLAB Drive/lib/casadi-folder'); 
import casadi.*


%%
x = MX.sym('x', 6); % [x y z gamma zeta v]
u = MX.sym('u', 3); % [omega_g omega_z a]
f_x_u = [x(6)*cos(x(4))*cos(x(5));
         x(6)*cos(x(4))*sin(x(5));
         x(6)*sin(x(4));
         u(1);
         u(2);
         u(3)];

X_dot = Function('X_dot',{x, u},{f_x_u});            % X_dot as function of X and U
clear x u f_x_u;
