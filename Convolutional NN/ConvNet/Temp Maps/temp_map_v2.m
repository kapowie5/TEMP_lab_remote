%% to create a circle (page 103 of Heat conduction, Latif M. Jiji, solution 108)
clear all
close all
clc

r = 0:1:100;  % half width of an image for T_full
r0 = r(end);
nr = length(r);
z = 0:1:120;    %Length of image
L = z(end);
T0 = 300;
Ta = 320;
q = 0.0000001;  %heat generation term (want near 0.01 or 0.1)
k_prop = 10;

for i_z = 1:length(z)
    for i_r = 1:length(r)  
        sum = 0;
        for k=1:100 %number of sumation iterations
            lambda(k) = (2*k-1)*pi/ (2*L);
            a(k) = (-2/(lambda(k)*L*besseli(0,(lambda(k)*r0))))*(T0-Ta+q/(k_prop*lambda(k)^2));
            term2(k) = besseli(0,(lambda(k)*r(i_r)))*sin(lambda(k)*z(i_z));
            sum = sum + a(k)*term2(k);
        end
        T(i_r,i_z) = T0 + (q*L^2/(2*k_prop))*(2*(z(i_z)/L) - (z(i_z)/L)^2) + sum;
    end
end

%%
figure;
plot(r,T(:,5))

figure;
plot(z,T(1,:))

figure;
 contourf(z,r,T)
 colorbar 
 
 %%
 clear r_full
T_full = [flipud(T);T(2:end,:)];
r_full = [fliplr(r),r(2:end)*-1];



figure;
 contourf(z,r_full,T_full)
 colorbar 
%%
% d_theta = 360/nr;
% theta  = [1:d_theta:360];

% [x,y] = pol2cart(theta,r);