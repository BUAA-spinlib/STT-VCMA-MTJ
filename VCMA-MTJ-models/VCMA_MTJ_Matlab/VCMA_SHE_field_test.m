
clc;
clear ;
sim_mid1 = 2000;   %Voltage start
sim_mid2 = 3000;   %Voltage stop, SHE start
sim_mid3 = 4000;   %SHE stop
sim_end = 10000;    %Simulation time
t_step = 1e-12;    %Simulation step
t_sim = 0:t_step:t_step*sim_end;

%------------------------Constants----------------------------%
uB = 9.274e-24;             %Bohr mageton in J/T
h_bar = 1.054e-34;          %Reduced Planck constant in Js
u0 = pi*4e-7;               %Vacuum permeability in H/m
e = 1.6e-19;                %Elementary charge in C
kb = 1.38e-23;              %Boltzman constant in J/K
gamma = 2*u0*uB/h_bar;      %Gyromagnetic ratio in m/As

%------------------------MTJ Parameters----------------------------%
a = 80e-9;                      %Width of the MTJ in m
b = 80e-9;                      %Length of the MTJ in m
tf = 1.1e-9;                    %Thickness of the freelayer in m                           
tox = 1.5e-9;                     %Thickness of the MgO barrier in m
P = 0.6;                        %Spin polarization of the tunnel current
alpha = 0.1;                   %Gilbert damping damping factor
T = 300;                        %Temperature in K
Ki = 1.0056364e-3;              %The anisotropy energy in J/m2
Ms = 1.2e6;                     %Saturation magnetization in A/m
ksi = 75e-15;                   %VCMA coefficient in J/vm
gammap = gamma/(1+alpha*alpha); %Reduced gyromagnetic ratio
v = tf*pi*b*a/4;                %Volume of free layer of the elliptical MTJ
delta = 40;                     %Thermal stability factor at V=0
TMR = 1;                        %TMR ratio at V=0,100%  
Rp = 5e4;                       %Magenetoresistance at parallel state, 100kOhm
A = (a^2)*pi/4;                 %MTJ area
eta = 0.3;                      %Spin hall angle
w = 100e-9;                    
l = 100e-9;
d = 3e-9;                      %Width,length and thichness of beta-W strip (heavy metal layer)
A2 = d*w;                     %Cross-sectional area of heavy metal layer 
rho = 200e-8;                   %Resistivity of beta-W
R2 = rho*l/(w*d);            %Resistance of beta-W
I_SHE = 60e-6;

%External filed;
Hx = 0;        
Hy = -8000;                     % An assisted field in y axis will speed up the switching, 100Oe
Hz = 0;   

%STT factor;
F = (gamma*h_bar)/(2*u0*e*tf*Ms);     

%Demagnetization filed;
Nx = 0.010613177892974;
Ny = 0.010613177892974;
Nz = 0.978773644214052;

%------------------------Initialization----------------------------%
phi = zeros(1,sim_end+1);
theta = zeros(1,sim_end+1);
power = zeros(1,sim_end+1);
phi(1) = 0;                              %initial azimuthal angle
theta(1) = pi-sqrt(1/(2*delta));         %initial polar angle
power(1) = 0;

for loop = 1:sim_mid1                               %Initial stage
    V = 0;                                          
    J_SHE = 0;
    Hk = (2*Ki)/(tf*Ms*u0)-(2*ksi*V)/(u0*Ms*tox*tf);  %The effective anisotropy field with VCMA effect  
    Ax = Hx-Nx*Ms*sin(theta(loop))*cos(phi(loop));        %x axis component of the effective magneitc field
    Ay = Hy-Ny*Ms*sin(theta(loop))*sin(phi(loop));        %y axis component of the effective magneitc field
    Az = Hz-Nz*Ms*cos(theta(loop))+Hk*cos(theta(loop));   %z axis component of the effective magneitc field
    dphi = gammap*(Ax*(-cos(theta(loop))*cos(phi(loop))-alpha*sin(phi(loop)))+Ay*(-cos(theta(loop))*sin(phi(loop))+alpha*cos(phi(loop)))+Az*sin(theta(loop)))/(sin(theta(loop)))+J_SHE*F*eta*(sin(phi(loop))-alpha*cos(phi(loop))*cos(theta(loop)))/(sin(theta(loop))*(1+alpha*alpha));           
    dtheta = gammap*(Ax*(alpha*cos(theta(loop))*cos(phi(loop))-sin(phi(loop)))+Ay*(alpha*cos(theta(loop))*sin(phi(loop))+cos(phi(loop)))-Az*alpha*sin(theta(loop)))-J_SHE*F*eta*(cos(phi(loop))*cos(theta(loop))+alpha*sin(phi(loop)))/(1+alpha*alpha);          
    R1 = Rp*(1+(V/0.5)^2+TMR)/(1+(V/0.5)^2+TMR*(1+cos(theta(loop)))/2);     %MTJ Resistance 
    dpower = V^2/R1+R2*(abs(J_SHE*A2))^2;                                  %Power consumption calculation
    phi(loop+1) = phi(loop)+t_step*dphi;                                    %Adopt the rectangular integral 
    theta(loop+1) = theta(loop)+t_step*dtheta;
    power(loop+1) = power(loop)+t_step*dpower;
end

for loop = sim_mid1:sim_mid2                        %VCMA and SHE
    V = 0.5;                                         %Applied votage for VCMA effect;
    J_SHE = -2.2e11;                                    %Applied current density for SHE
    Hk = (2*Ki)/(tf*Ms*u0)-(2*ksi*V)/(u0*Ms*tox*tf);  %The effective anisotropy field with VCMA effect  
    Ax = Hx-Nx*Ms*sin(theta(loop))*cos(phi(loop));        %x axis component of the effective magneitc field
    Ay = Hy-Ny*Ms*sin(theta(loop))*sin(phi(loop));        %y axis component of the effective magneitc field
    Az = Hz-Nz*Ms*cos(theta(loop))+Hk*cos(theta(loop));   %z axis component of the effective magneitc field
    dphi = gammap*(Ax*(-cos(theta(loop))*cos(phi(loop))-alpha*sin(phi(loop)))+Ay*(-cos(theta(loop))*sin(phi(loop))+alpha*cos(phi(loop)))+Az*sin(theta(loop)))/(sin(theta(loop)))+J_SHE*F*eta*(sin(phi(loop))-alpha*cos(phi(loop))*cos(theta(loop)))/(sin(theta(loop))*(1+alpha*alpha));           
    dtheta = gammap*(Ax*(alpha*cos(theta(loop))*cos(phi(loop))-sin(phi(loop)))+Ay*(alpha*cos(theta(loop))*sin(phi(loop))+cos(phi(loop)))-Az*alpha*sin(theta(loop)))-J_SHE*F*eta*(cos(phi(loop))*cos(theta(loop))+alpha*sin(phi(loop)))/(1+alpha*alpha);          
    R1 = Rp*(1+(V/0.5)^2+TMR)/(1+(V/0.5)^2+TMR*(1+cos(theta(loop)))/2);     %MTJ Resistance 
    dpower = V^2/R1+R2*(abs(J_SHE*A2))^2;                                  %Power consumption calculation
    phi(loop+1) = phi(loop)+t_step*dphi;                                    %Adopt the rectangular integral 
    theta(loop+1) = theta(loop)+t_step*dtheta;
    power(loop+1) = power(loop)+t_step*dpower;
end

for loop = sim_mid2:sim_mid3                        %SHE
    V = 0;                                          
    J_SHE = 0;
    Hk = (2*Ki)/(tf*Ms*u0)-(2*ksi*V)/(u0*Ms*tox*tf);  %The effective anisotropy field with VCMA effect  
    Ax = Hx-Nx*Ms*sin(theta(loop))*cos(phi(loop));        %x axis component of the effective magneitc field
    Ay = Hy-Ny*Ms*sin(theta(loop))*sin(phi(loop));        %y axis component of the effective magneitc field
    Az = Hz-Nz*Ms*cos(theta(loop))+Hk*cos(theta(loop));   %z axis component of the effective magneitc field
    dphi = gammap*(Ax*(-cos(theta(loop))*cos(phi(loop))-alpha*sin(phi(loop)))+Ay*(-cos(theta(loop))*sin(phi(loop))+alpha*cos(phi(loop)))+Az*sin(theta(loop)))/(sin(theta(loop)))+J_SHE*F*eta*(sin(phi(loop))-alpha*cos(phi(loop))*cos(theta(loop)))/(sin(theta(loop))*(1+alpha*alpha));           
    dtheta = gammap*(Ax*(alpha*cos(theta(loop))*cos(phi(loop))-sin(phi(loop)))+Ay*(alpha*cos(theta(loop))*sin(phi(loop))+cos(phi(loop)))-Az*alpha*sin(theta(loop)))-J_SHE*F*eta*(cos(phi(loop))*cos(theta(loop))+alpha*sin(phi(loop)))/(1+alpha*alpha);          
    R1 = Rp*(1+(V/0.5)^2+TMR)/(1+(V/0.5)^2+TMR*(1+cos(theta(loop)))/2);     %MTJ Resistance 
    dpower = V^2/R1+R2*(abs(J_SHE*A2))^2;                                  %Power consumption calculation
    phi(loop+1) = phi(loop)+t_step*dphi;                                    %Adopt the rectangular integral 
    theta(loop+1) = theta(loop)+t_step*dtheta;
    power(loop+1) = power(loop)+t_step*dpower;
end

for loop = sim_mid3:sim_end                        %Stabilization stage
    V = 0;                                         
    J_SHE = 0;
    Hk = (2*Ki)/(tf*Ms*u0)-(2*ksi*V)/(u0*Ms*tox*tf);  %The effective anisotropy field with VCMA effect  
    Ax = Hx-Nx*Ms*sin(theta(loop))*cos(phi(loop));        %x axis component of the effective magneitc field
    Ay = Hy-Ny*Ms*sin(theta(loop))*sin(phi(loop));        %y axis component of the effective magneitc field
    Az = Hz-Nz*Ms*cos(theta(loop))+Hk*cos(theta(loop));   %z axis component of the effective magneitc field
    dphi = gammap*(Ax*(-cos(theta(loop))*cos(phi(loop))-alpha*sin(phi(loop)))+Ay*(-cos(theta(loop))*sin(phi(loop))+alpha*cos(phi(loop)))+Az*sin(theta(loop)))/(sin(theta(loop)))+J_SHE*F*eta*(sin(phi(loop))-alpha*cos(phi(loop))*cos(theta(loop)))/(sin(theta(loop))*(1+alpha*alpha));           
    dtheta = gammap*(Ax*(alpha*cos(theta(loop))*cos(phi(loop))-sin(phi(loop)))+Ay*(alpha*cos(theta(loop))*sin(phi(loop))+cos(phi(loop)))-Az*alpha*sin(theta(loop)))-J_SHE*F*eta*(cos(phi(loop))*cos(theta(loop))+alpha*sin(phi(loop)))/(1+alpha*alpha);          
    R1 = Rp*(1+(V/0.5)^2+TMR)/(1+(V/0.5)^2+TMR*(1+cos(theta(loop)))/2);     %MTJ Resistance 
    dpower = V^2/R1+R2*(abs(J_SHE*A2))^2;                                  %Power consumption calculation
    phi(loop+1) = phi(loop)+t_step*dphi;                                    %Adopt the rectangular integral 
    theta(loop+1) = theta(loop)+t_step*dtheta;
    power(loop+1) = power(loop)+t_step*dpower;
end
phi_angle = mod(phi,2*pi)*180/pi;               
theta_angle = mod(theta,2*pi)*180/pi;

mx = sin(theta).*cos(phi);              %X axis component of magnetization vector m of free layer
my = sin(theta).*sin(phi);              %Y axis component of magnetization vector m of free layer
mz = cos(theta);                        %Z axis component of magnetization vector m of free layer

figure(1)
plot(t_sim,phi_angle,'g');
title('\phi');

figure(2)
plot(t_sim,theta_angle,'b');
title('\theta');

figure(3)
%{
plot(t_sim,mx,'k','linewidth',1.0);
title('mx');
hold on;
plot(t_sim,my,'b','linewidth',1.0);
title('my');
%}
hold on;
plot(t_sim,mz,'r','linewidth',2.0);
title('mz');
hold on;

figure(4)
plot(t_sim,power,'r','linewidth',2.0);
title('power');
hold on;

figure(5)
plot3(my,mx,mz,'r','linewidth',2.0);
grid on; 
axis([-1 1 -1 1 -1 1]);

