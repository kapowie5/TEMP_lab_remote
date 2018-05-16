%THIS IS PETER
% https://www.mathworks.com/help/matlab/ref/integral.html
clc
clear all
close all

rho=2328; % kg/m^3                                   %Peter: 10960 to 2328 for silicon
cp=712; % J/kg*K                                      %237 to 712
k=150; % W/m*K                                          %8 to 150
alpha = k/(cp*rho);
kg = 0.024; % air W/m*K
alphag = 1.9e-5; % air m^2/s
x1 = 418.95e-6; % m
x2 = 438.88e-6; % m
L = 450e-6; % m cantilever length
D1 = 31e-6; % Spot size, m
D2 = 65.3e-6; % spot size, m
thick = 2.5e-6; % m
h = 22e-6;
vol = D1*D2*thick;
P0 = 7.9e-3;
g = P0/vol;
eps = 0.1;          % emissivity
sigma = 5.67e-8;    % Stephan-Boltzman Constant
T0 = 300;
f0 = 10^5;

%% Freq Resolved Raman
error_series=1e-15;
sum_series=0;
delta_T=1;
clear m
m=1;
clear t
t = [0:0.0001:0.03]; % Time
l = L-h/2;

while (delta_T>error_series)
        Cm = ((1-(-1)^m)*cos((m*pi*x1)/(2*l)))^2/(m^4*pi^4);                %got rid of a negative in the cos
        delta_T=Cm * ((1-exp(-(m^2*pi^2*alpha*t)/(4*l^2)))/...
            (1+exp(-(m^2*pi^2*alpha)/(8*f0*l^2))));
        sum_series=sum_series+delta_T;
        m=m+1;
end
C0 = 8*g*L^3/((l-x1)*k);
theta_pulse = C0 * sum_series;
T = theta_pulse+T0;

figure;
plot(t,T);
xlabel('Time (s)')
ylabel('Temp (K)')
title('Frequency Resolved Raman')