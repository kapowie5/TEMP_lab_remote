% function y = computeRaquel(delta,r,a,w,alpha,G,H)
%     beta = delta^2 + 1i*w/alpha;
%     y = exp(-delta^2*a^2/8)/beta * 1/(1+G+H) * besselj(0,delta*r) * delta;
% end 
% 
% f = @computRaquel;
% q = integral(f,0,inf,'RelTol',1e-8,'AbsTol',1e-13)  
a = 10e-6;
h = 6.3;
ks = 10;
G = 0;
r=15e-6;
alpha = 5e-7;
alphag = 1e-5;
w = 100;
kg = 0;

% kg*(d.^2 + 1i*w./alphag)./(k*(d.^2 + 1i*w./alpha));

fun = @(d)d.*besselj(0,d*r).*exp(-d.^2*a^2/8)./(d.^2 + 1i*w/alpha)* 1./(1+G+h./(ks*(d.^2 + 1i*w./alpha)));
fun = @(d)d.*besselj(0,d*r).*exp(-d.^2*a^2/8)./(d.^2 + 1i*w/alpha)*...
    1./(1 + (d.^2 + 1i*w./alphag)./(ks*(d.^2 + 1i*w./alpha)) + h./(ks*(d.^2 + 1i*w./alpha)));

% fun2 = @(d)d.*(d.^2 + 1i*w./alphag)./(ks*(d.^2 + 1i*w./alpha)).*kg;
% besselj(0,t)
q = integral(fun,0,Inf,'RelTol',1e-8,'AbsTol',1e-13)  
% q2 = integral(fun2,0,Inf,'RelTol',1e-8,'AbsTol',1e-13)  

% close all
% y=fun(0:1);
% plot(0:1,abs(y))
