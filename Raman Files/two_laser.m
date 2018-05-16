clc
%% Two laser model - w/o rad
% http://dx.doi.org/10.1063/1.4867166 View
figure; hold on
r=logspace(-7,-5);
T0 = 300; % K
%P_abs = 10e-3; % W
d = 250e-9;
for count = 1:7
    if count==1
        P_abs = 2.7e-3;
        P_abs = P_abs*.51455;  %http://www2.ensc.sfu.ca/~glennc/e894/e894l15g.pdf A=1-R
    end
    if count==2
        P_abs = 4.7e-3;
        P_abs = P_abs*.51455;
    end
    if count==3
        P_abs = 6e-3;
        P_abs = P_abs*.51455;
    end
    if count==4
        P_abs = 11.4e-3;
        P_abs = P_abs*.51455;
    end
    if count==5
        P_abs = 16.8e-3;
        P_abs = P_abs*.51455;
    end
    if count==6
        P_abs = 26.8e-3;
        P_abs = P_abs*.51455;
    end
    if count==7
        P_abs = 32.2e-3;
        P_abs = P_abs*.51455;
    end
    t = .25e-6; % m               thikness, (d in literature)
    r0 = r(end);
    k=80; % W/m*K    for silicon

    T_noRad = T0 - P_abs/(2*pi*t*k)*log(r/r0);
%     figure;
    plot(r,T_noRad)
    
%     title('Two Laser w/o rad Raman')
%     xlabel('r (m)')
%     ylabel('Temp (K)')
end
% t = .25e-6; % m               thikness, (d in literature)
% r0 = r(end);
% k=150; % W/m*K    for silicon
% T_noRad = T0 - P_abs/(2*pi*t*k)*log(r/r0);
% figure;plot(r,T_noRad)
title('Two Laser w/o rad Raman')
xlabel('r (m)')
ylabel('Temp (K)')


% figure;semilogx(r,T_noRad)
% figure;semilogy(r,T_noRad)

%%
%% Two laser model - w/ rad
% close all
% From Incropera C .17

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
    c1(i) = (-0.5*q*r1/k)*(besselk (0,N*r2)/(N*besseli(1,N*r1)*besselk(0,N*r2) +...
        N*besseli(0,N*r2)*besselk(1,N*r1)));
    c2(i) = -c1(i) * besseli (0,N*r2)/besselk(0,N*r2);
    T2(i) = c1(i)*besseli(0,N*r(i)) + c2(i)*besselk(0,N*r(i)) + T0;
    T2_r1 = c1(i)*besseli(0,N*r1) + c2(i)*besselk(0,N*r1) + T0;;
    T1(i) = ((q*r1^2)/(4*k))*(1-r (i)^2/r1^2) + T2_r1;
    
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
%     T1(i) = ((q*r1^2)/(4*k))*(1-r (i)^2/r1^2) + T2_r1;
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
% % From Incropera C .17
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
% A = sqrt(z)*(sinh (sqrt(z)*log10(r2))^2 - cosh(sqrt(z)*log10(r2))*cosh(sqrt(z)*log10(r1)))/sinh(sqrt(z)*log10(r2));
% 
% for i =1:length(r)
%     B(i) = ((q*r1^2)/(4*k))*(1-r (i)^2/r1^2);
%     C(i) = ((q*r1^2)/(2*k))*A^(-1)*(cosh(sqrt(z)*log10(r1))+cosh(sqrt(z)*log10(r2))*cosh (sqrt(z)*log10(r1))/sinh(sqrt(z)*log10(r2)));
%     D(i) = ((q*r1^2)/(2*k))*A^(-1)*cosh(sqrt(z)*log10(r(i)));
%     E(i) = ((1*r1^2)/(2*k))*A^(-1)*cosh (sqrt(z)*log10(r2))/sinh(sqrt(z)*log10(r2))*sinh(sqrt(z)*log10(r(i)));
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