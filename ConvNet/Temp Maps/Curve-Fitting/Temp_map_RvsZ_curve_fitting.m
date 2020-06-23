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
        T(i_r,i_z) = T0 + (q*L^2/(2*k_prop))*(2*(z(i_z)/L) - (z(i_z)/L)^2) + sum + normrnd(0,0.1); 
    end
end

%%
figure;
plot(r,T(:,5))

figure;
plot(z,T(1,:))

figure;
 contourf(z,r,T)
 
 %%
 clear r_full
 fsize = 16;
True_Temp = [flipud(T);T(2:end,:)];
r_full = [fliplr(r),r(2:end)*-1];



figure;
 contourf(z,r_full,True_Temp) 
 h = colorbar; set(get(h,'label'),'string','Temp (K)','FontWeight','bold','FontSize',fsize);
 set(gca,'FontWeight','bold','FontSize',fsize,'XTickLabel',{},'YTickLabel',{}); 
 saveas(gcf,'Temp Map - Temps (RvsZ).fig')
print('Temp Map - Temps (RvsZ)','-dpng','-r300')
%%
% d_theta = 360/nr;
% theta  = [1:d_theta:360];

% [x,y] = pol2cart(theta,r);

%% 
% This creates a linear relationship between the temperature and pixel
% values of an image, and assumes a 16 bit pixel depth so that smaller
% temperature variations can be seen (a lot of greyscale cameras are 12
% bit, but you can change the bit depther easily
% look at Linear to Curve fit to NN comparison (Autosaved).xlsx


%Image 16 bit grayscale intensity due to temperature (290-380 range)
m =     [-8.24	-54.03	-210.27	-272.26	-66.21 ;
        3067.50	26378.00	105405.26	139946.00	36360.00];
        
% %Image 16 bit grayscale intensity due to temperature (300-312 range)
% m =     [-61.78	-405.22	-1577.00	-2041.90	-496.58 ;
%         19213.00	132276.00	517539.07	673575.00	166133.00];

Izero   = uint16(zeros([size(True_Temp),3]));
Iblue   = uint16(m(1,1)*True_Temp + m(2,1));
Igreen  = uint16(m(1,2)*True_Temp + m(2,2));
Iyellow = uint16(m(1,3)*True_Temp + m(2,3));
Iorange = uint16(m(1,4)*True_Temp + m(2,4));
Ired    = uint16(m(1,5)*True_Temp + m(2,5));

%%

Red_vec=Ired(:);
Orange_vec=Iorange(:);
Yellow_vec=Iyellow(:);
Green_vec=Igreen(:);
Blue_vec=Iblue(:);
T_vec=True_Temp(:);

%%
figure;
plot(T_vec,Red_vec,'r',T_vec,Orange_vec,'k',T_vec,Yellow_vec,'y',T_vec,Green_vec,'g',T_vec,Blue_vec,'b')



% %%
% % close all
% imBlue = Izero; imBlue(:,:,3) = Iblue*20;
% imGreen = Izero; imGreen(:,:,2) = Igreen;
% imYellow = Izero; imYellow(:,:,1) = 204*256;       imYellow(:,:,2) = 204*256;     imYellow(:,:,3) = Iyellow/1.5;
% imOrange = Izero; imOrange(:,:,1) = 2^16;       imOrange(:,:,2) = Iorange/2;
% imRed = Izero; imRed(:,:,1) = Ired/1.5;
% 
% figure;imshow(imBlue)
% figure;imshow(imGreen)
% figure;imshow(imYellow)
% figure;imshow(imOrange)
% figure;imshow(imRed)
% 
% % imViolet = Izero; imOrange(:,:,1) = Iorange;    imOrange(:,:,3) = Iorange;
% % figure;imshow(imViolet)