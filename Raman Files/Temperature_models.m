% https://www.mathworks.com/help/matlab/ref/integral.html
clc
clear all
close all

rho=10960; %kg/m^3
cp=237; %J/kg*K
k=8; %W/m*K
alpha = k/(cp*rho);
kg = 0.024; %air W/m*K
alphag = 1.9e-5; %air m^2/s
x1 = 418.95e-6; %m
x2 = 438.88e-6; %m
L = 450e-6; %m cantilever length
D1 = 31e-6; %Spot size, m
D2 = 65.3e-6; %spot size, m
thick = 2.5e-6; %m
h = 22e-6;
vol = D1*D2*thick;
P0 = 7.9e-3;
g = P0/vol;
eps = 0.1;          %emissivity
sigma = 5.67e-8;    %Stephan-Boltzman Constant
T0 = 300;
f0 = 2;

%% Time domain Raman
error_series=1e-15;
sum_series=0;
delta_T=1;
m=1;
t = [0:0.1:10]; %Time

while (delta_T>error_series)
        delta_T=(1/m^4/pi^4)*(1-exp(-m^2*pi^2*alpha*t))*...
            (cos(m*pi*x1/L) - cos(m*pi*x2/L))^2;
        sum_series=sum_series+delta_T;
        m=m+1;
    end

% top = (1/m^4/pi^4)*(1-exp(-m^2*pi^2*alpha*t))*(cos(m*pi*x1/L) - cos(m*pi*x2/L))^2;
% bottom = (1/m^4/pi^4)*(cos(m*pi*x1/L) - cos(m*pi*x2/L))^2;
% theta_s = top/bottom;

% sum = (1/m^4/pi^4)*(1-exp(-m^2*pi^2*alpha*t))*(cos(m*pi*x1/L) - cos(m*pi*x2/L))^2;
theta = (2*g*L^3)/((x2-x1)*k) * sum_series;
T = theta + T0;

figure;
plot(t,T);
xlabel('Time (s)')
ylabel('Temp (K)')
title('Time Domain Raman')


%% Freq Resolved Raman
error_series=1e-15;
sum_series=0;
delta_T=1;
m=1;
clear t
t = [0:0.001:0.1]; %Time
l = L-h/2;

while (delta_T>error_series)
        Cm = ((1-(-1)^m)*cos(-(m*pi*x1)/(2*l)))^2/(m^4*pi^4);
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
Fo = alpha*t/L;
error_series=1e-15;
sum_series=0;
delta_T=1;
mm=1;
q = g;

Tss = (1/L^2/He^2)*(1-(sinh(L*He*xs)+sinh(L*He*(1-xs)))/sinh(L*He));
while (delta_T>error_series)
        delta_T=exp(-((2*mm-1)^2*pi^2 + L^2*He^2) * Fo) * sin((2*m-1)*pi*xs) / ...
            ((2*m-1)^3*pi^3 + L^2*He^2*(2*m-1)*pi);
        sum_series=sum_series+delta_T;
        mm=mm+1;
    end
Tts = -4*sum_series;
Ts = Tss + Tts;
T = Ts*q*L^2/k + T0;


figure;
plot(t,T);
xlabel('Time (s)')
ylabel('Temp (K)')
title('Spider Silk Raman')

%% Raquel - Thin Film, same surface temperature
P0 = 1e-3; %W of laser power
a = 10e-6; %m, diameter laser
h = 6.3;
ks = 10;
G = 0;
r=15e-6;
alpha = 5e-7;
% w = [1:100];
omega = [1:100];
% w = 1;

clear Ts
for ii=1:100
    w = omega(ii);

fun = @(d)d.*besselj(0,d*r).*exp(-d.^2*a^2/8)./(d.^2 + 1i*w/alpha)*...
    1./(1 + (d.^2 + 1i*w./alphag)./(ks*(d.^2 + 1i*w./alpha)) + h./(ks*(d.^2 + 1i*w./alpha)));
q = integral(fun,0,Inf,'RelTol',1e-8,'AbsTol',1e-13);
Ts(ii) = P0/(4*pi*k) * q;
end

figure;
plot(omega,abs(Ts));
xlabel('\omega (rad/s)')
ylabel('Temp Amplitude')
title('Raquel Film Magnitude')
figure;
plot(omega,unwrap(angle((Ts))));
xlabel('\omega (rad/s)')
ylabel('Temp Phase (rad)')
title('Raquel Film Phase')


%% Zilong - Thin Film







































%% Two laser model - w/o rad
% http://dx.doi.org/10.1063/1.4867166 View

r=logspace(-7,-5);
T0 = 300; %K
P_abs = 10e-3; %W
t = 3e-6; %m
r0 = r(end);

T_noRad = T0 - P_abs/(2*pi*t*k)*log(r/r0);
figure;plot(r,T_noRad)
title('Two Laser w/o rad Raman')
xlabel('r (m)')
ylabel('Temp (K)')
% figure;semilogx(r,T_noRad)
% figure;semilogy(r,T_noRad)

%%
%% Two laser model - w/ rad
% close all
% From Incropera C.17

r2 = r(end);
r1 = 10e-6;
epi = 0.3;
sigma = 5.67e-8;
h_r = 4*epi*sigma*T0^3;
% h_r=0;
q = P_abs/(pi*t*r1);
N = sqrt(2*h_r/(k*t));
I = find(r < r1);
Ir = I(end);

%%

for i =1:length(r)
    c1(i) = (-0.5*q*r1/k)*(besselk(0,N*r2)/(N*besseli(1,N*r1)*besselk(0,N*r2) +...
        N*besseli(0,N*r2)*besselk(1,N*r1)));
    c2(i) = -c1(i) * besseli(0,N*r2)/besselk(0,N*r2);
    T2(i) = c1(i)*besseli(0,N*r(i)) + c2(i)*besselk(0,N*r(i)) + T0;
    T2_r1 = c1(i)*besseli(0,N*r1) + c2(i)*besselk(0,N*r1) + T0;;
    T1(i) = ((q*r1^2)/(4*k))*(1-r(i)^2/r1^2) + T2_r1;
    
    if (i<=Ir)
        T_Rad(i) = T1(i);
    else
        T_Rad(i) = T2(i);
    end
end
% for i =1:length(r)
%     A(i) = -q*r1/(2*k*N*besselk(1,N*r1));
%     B(i) = ((besseli(0,N*r2)*besselk(1,N*r1) + besseli(1,N*r1)*besselk(0,N*r2))/besselk(1,N*r1))^-1;
%     c1(i) = A(i)*B(i);
%     C(i) = (q*r1 + 2*k*c1(i)*N*besseli(1,N*r1))/(2*k);
%     T2(i) = A(i)*B(i)*besseli(0,N*r(i)) + C(i)*besselk(0,N*r(i)) + T0;
%     T2_r1 = A(i)*B(i)*besseli(0,N*r1) + C(i)*besselk(0,N*r1) + T0;
%     T1(i) = ((q*r1^2)/(4*k))*(1-r(i)^2/r1^2) + T2_r1;
%     
%     if (i<=Ir)
%         T_Rad(i) = T1(i);
%     else
%         T_Rad(i) = T2(i);
%     end
% end

% 
% figure;loglog(r,T_Rad)
% figure;semilogx(r,T_Rad)
% figure;semilogy(r,T_Rad)

% figure;loglog(r,T_Rad)
% figure;semilogx(r,T_Rad)
% figure;semilogy(r,T_Rad)

% %% Two laser model - w/ rad
% % close all
% % From Incropera C.17
% 
% r2 = r(end);
% r1 = 1e-6;
% epi = 0.1;
% sigma = 5.67e-8;
% h_r = 4*epi*sigma*T0^3;
% % h_r=0;
% q = P_abs/(2*pi*t);
% z = sqrt(2*h_r/(k*t));
% 
% A = sqrt(z)*(sinh(sqrt(z)*log10(r2))^2 - cosh(sqrt(z)*log10(r2))*cosh(sqrt(z)*log10(r1)))/sinh(sqrt(z)*log10(r2));
% 
% for i =1:length(r)
%     B(i) = ((q*r1^2)/(4*k))*(1-r(i)^2/r1^2);
%     C(i) = ((q*r1^2)/(2*k))*A^(-1)*(cosh(sqrt(z)*log10(r1))+cosh(sqrt(z)*log10(r2))*cosh(sqrt(z)*log10(r1))/sinh(sqrt(z)*log10(r2)));
%     D(i) = ((q*r1^2)/(2*k))*A^(-1)*cosh(sqrt(z)*log10(r(i)));
%     E(i) = ((1*r1^2)/(2*k))*A^(-1)*cosh(sqrt(z)*log10(r2))/sinh(sqrt(z)*log10(r2))*sinh(sqrt(z)*log10(r(i)));
%     T_Rad(i) = B(i)+C(i)+D(i)+E(i) + T0;
% end
% 
% 
% % figure;loglog(r,T_Rad)
% % figure;semilogx(r,T_Rad)
% % figure;semilogy(r,T_Rad)

%%
% figure;
% semilogx(r,T_Rad)
% hold all
% semilogx(r,T_noRad,'+')
% legend('Rad','No Rad')

% 
% figure;
% loglog(r,T_Rad)
% hold all
% loglog(r,T_noRad,'+')
% legend('Rad','No Rad')