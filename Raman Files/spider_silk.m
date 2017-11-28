clc
clear all
close all

%rho=10960; % kg/m^3
%cp=237; % J/kg*K
k=1.25; % W/m*K                  %Peter: values from "Termal characterization of natural and synthetic spider silks (Munro)
alpha = 5.7e-7;%k/(cp*rho);        %Peter: values interpolated due to thickness and diameter used in paper
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
f0 = 2;

%% Spider Silk Model
% A_c = D^2*pi/4;
A_c = D1 * thick;
% A_l = 2*pi*D*L;
A_l = L*thick;
hr = 4*T0*eps*sigma;
hc = 0;
He = sqrt((hr+hc)*A_l/(k*A_c*L));
x = 100e-6;
xs = x/L;
t = [0:0.1:10]; % Time                     
Fo = alpha*t/L;
error_series=1e-15;
sum_series=0;
delta_T=1;
mm=1;
q = g;

Tss = (1/L^2/He^2)*(1-(sinh(L*He*xs)+sinh(L*He*(1-xs)))/sinh(L*He));
while (delta_T>error_series)
        delta_T=exp(-((2*mm-1)^2*pi^2 + L^2*He^2) * Fo) * sin ((2*mm-1)*pi*xs) / ((2*mm-1)^3*pi^3 + L^2*He^2*(2*mm-1)*pi);
        sum_series=sum_series+delta_T;
        mm=mm+1;                          %Peter: the mm in the loop was m in some locations
end
Tts = -4*sum_series;
Ts = Tss + Tts;
T = Ts*q*L^2/k + T0;

figure;
plot(t,T);
xlabel('Time (s)')
ylabel('Temp (K)')
title('Spider Silk Raman')