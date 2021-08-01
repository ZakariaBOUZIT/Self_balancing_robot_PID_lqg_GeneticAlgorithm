%credit to Steve Brunton 
close all; clear all; clc
M = 1;    % mass of the chassis  
m = 0.2;  % mass of the wheels and shaft
b = 0.1;  % estimate of viscous friction coefficient (N-m-s)
I = 0.0005;% moment of inertia of the pendulum
g = 9.8;  %acceleration due to gravity (m/s^2)
l = 0.125;  %length to pendulum center of mass
p = I*(M+m)+M*m*l^2; %denominator for the A and B matrices
dt = 0.001;
PopSize = 25;
MaxGenerations = 10;
s = tf('s');
G = s+1/s^2+1; 
options =  optimoptions(@ga,'PopulationSize',PopSize,'MaxGenerations',MaxGenerations);
[x,fval] = ga(@(K)pidtest(G,dt,K),3,-eye(3),zeros(3,1))
function J = pidtest(G,dt,parms)
    s = tf('s');
    K = parms(1) + parms(2)/s + parms(3)*s/(1+.001*s);
    Loop = series(K,G);
    ClosedLoop = feedback(Loop,1);
    t = 0:dt:20;
    [y,t] = step(ClosedLoop,t);
    CTRLtf = K/(1+K*G);
    u = lsim(CTRLtf,1-y,t);
    Q = 1;
    R = .001;
    J = dt*sum(Q*(1-y(:)).^2+R*u(:).^2)
    [y,t] = step(ClosedLoop,t);
    plot(t,y,'LineWidth',2,'color','r')
    drawnow
end
