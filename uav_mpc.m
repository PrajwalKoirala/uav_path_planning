clear all; close all;

addpath('/MATLAB Drive/lib/casadi-folder'); 
import casadi.*


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
%%
% Variables and parameters
N = 20; dt = 0.5;
simtime = 15;
opti=casadi.Opti();
% X = opti.variable(6 ,N+1);
% U = opti.variable(3, N);
x = opti.variable(1, N+1);      y = opti.variable(1, N+1);     z = opti.variable(1, N+1);
gamma = opti.variable(1, N+1);  zeta = opti.variable(1, N+1);  v = opti.variable(1, N+1);
X = [x; y; z; gamma; zeta; v];
omega_g = opti.variable(1, N);
omega_z = opti.variable(1, N);
acc = opti.variable(1, N);
U = [omega_g; omega_z; acc];
X_0 = opti.parameter(6,1);
X_f = opti.parameter(6,1);

%Obstacle:
No_obs = 2; 
Obs_pos = opti.parameter(2,No_obs); Obs_rad = opti.parameter(1,No_obs)+1;
obst_cost = 0;
for i = 1:No_obs
    h =  log( ((x-Obs_pos(1,i))/Obs_rad(i)).^2 + ((y-Obs_pos(2,i))/Obs_rad(i)).^2  );
    opti.subject_to(h>0);  
    obst_cost = obst_cost + sum(sum(exp(7*exp(-h))));
end

%%
opti.minimize( 0.02*sum(sum((x-X_f(1)).^2)) + 0.3*sum(sum((y).^2)) +  0.2*sum(sum((z-X_f(3)).^2)) + obst_cost );  % 0.02*sum(sum((x-X_f(1)).^2)) + 0.5*sum(sum((x-y).^2)) +  0.5*sum(sum((z-X_f(3)).^2)) 
for k = 1:N
    opti.subject_to( X(:,k+1)== X(:,k) + dt*X_dot(X(:,k),U(:,k)) );
end

opti.subject_to(X(:,1) == X_0); % initial conditions 
opti.subject_to(v <= 30); opti.subject_to(15 <= v);  % velocity constraints
opti.subject_to(acc <= 3); opti.subject_to(-3 <= acc);  % acceleration constraints
opti.subject_to(gamma <= pi/4); opti.subject_to(-pi/4 <= gamma);  % pitch angle constraints
opti.subject_to(zeta <= pi); opti.subject_to(-pi <= zeta);  % yaw angle constraints
opti.subject_to(z <= 100); opti.subject_to(0 <= z);  % height constraints
opti.subject_to(omega_z <= 0.5); opti.subject_to(-0.5 <= omega_z);
opti.subject_to(omega_g <= 0.5); opti.subject_to(-0.5 <= omega_g);

%
Obs_rad = Obs_rad+1;
opti.subject_to(h>0);  

p_opts = struct('expand',true);
s_opts = struct('max_iter',1000);
opti.solver('ipopt',p_opts, s_opts);

%%
% x_to_u = opti.to_function('x_to_u',{X_0, X_f}, {U(:,1)});

%% 
X_now = [0;0;0;0;0;15];
X_ref = [600;0;40;0;0;0];
iter = simtime/dt;
X_es = [X_now, zeros(6,iter)];
U_s = zeros(3,iter); 
% opti.set_initial(U,repmat(rand(3,1),1,N)); 
% opti.set_initial(X, X_now+linspace(0,1,N+1).*(X_ref-X_now)); 
%
Obst(:,1) = [200;5;10];    % [x_c;y_c;rad]
Obst(:,2) = [260;-10;10];
opti.set_value(Obs_pos, Obst(1:2,:)); opti.set_value(Obs_rad, Obst(3,:));

%
tic
for i = 1:iter
    opti.set_value(X_0, X_now); opti.set_value(X_f, X_ref);
    sol = opti.solve();
    u = sol.value(U); x_set = sol.value(X);
    opti.set_initial(U,u); opti.set_initial(X,x_set);
    u = u(:,1);
    U_s(:,i) = u;  
    X_now = X_now + dt* full(X_dot(X_now, u));
    X_es(:,i+1) = X_now;
end
toc
%%
close all;
x = X_es(1,:);
y = X_es(2,:);
z = X_es(3,:);
gamma = X_es(4,:);
zeta = X_es(5,:);
v = X_es(6,:);

omega_g = U_s(1,:);
omega_z = U_s(2,:);
acc = U_s(3,:);
t = [0:iter]*dt;
figure; hold on;
plot(t, x)
plot(t, y)
plot(t, z)
legend({'x','y','z'})
figure; hold on;
plot(t, v)
stairs(t, [acc, nan])
legend({'velocity','acceleration'})

figure; hold on;
stairs(t, [omega_g, nan])
stairs(t, [omega_z, nan])
legend({'omega_g','omega_z'})

figure; axis equal; hold on;
plot(x, y);  %(x,y)
for i = 1: No_obs
    theta = linspace(0,2*pi,100); plot(Obst(1,i) + Obst(3,i)*cos(theta), Obst(2,i) + Obst(3,i)*sin(theta));
end


figure; axis equal; hold on;
plot3(x, y, z);  view([60 25]);
for i = 1: No_obs
    [X,Y,Z] = cylinder(Obst(3,i));
    Z = Z*60-30;
    surf( X+Obst(1,i), Y+Obst(2,i), Z+X_ref(3) );
end


%%
% figure; hold on; axis([0 700 -5 5]);
% plot(x, y);
% for  i = 1:length(x)
%     plot(x(i),y(i),'o')
%     pause(0.02)
% end
