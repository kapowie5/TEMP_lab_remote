clc
%% Raquel - Film, same surface temperature
P0 = 1e-3; % W of laser power
a = 500e-6; % m, diameter laser
h = 50;
ks = 0.3;
kg = 0.026;
G = 0;
alpha = .1e-6;
alphag = 2.2e-5;
f = 0.05;
w = 2*pi*f; %figure 3.4
r = -0.0024:0.00001:0.0024;
clear Ts

for ii=1:length(r) 
    %equation 3.3; besselj matlab function J; 1i returns inaginary i;
    fun = @(d) d.*besselj (0,d*r(ii)).*exp (-(d.^2*a.^2)./8)./(d.^2 + 1i*w/alpha)*...
        1./(1 + (kg*(d.^2 + 1i*w./alphag)./(ks*(d.^2 + 1i*w./alpha))) + (h./(ks*(d.^2 + 1i*w./alpha))));
    q = integral(fun,0,Inf,'RelTol',1e-8,'AbsTol',1e-13); %relative and absolute tolerances
    Ts(ii) = P0/(4*pi*ks) * q;    
end

figure;
%plot(omega,abs(Ts));
plot(r,abs(Ts));
%plot(r,log(abs(Ts).*r),'c*');  %to match figure 3.4
%xlabel('\omega (rad/s)')
xlabel('radius (m)')
ylabel('Temp Amplitude')
title('Raquel Film Magnitude')
figure;
%plot(omega,unwrap(angle((Ts))));
plot(r,unwrap(angle((Ts)))-2*pi);
%xlabel('\omega (rad/s)')
xlabel('radius (m)')
ylabel('Temp Phase (rad)')
title('Raquel Film Phase')