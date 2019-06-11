clear all
close all
fsize = 16;

%Vary These temperatures to make different images
T1 = 290; %Top and bottom temperatures (K)
T2 = 340; %Side Temperatures (K)

%Dimensions of your image
% L = 121*4;  %pixels
% W = 201*4;  %pixels

L = 1080;  %pixels
W = 1920;  %pixels


x = 0:1:L-1;
y = 0:1:W-1;
maxnum = 1;
Igreen = zeros(length(y),length(x));

Iblue = zeros(length(y),length(x));

Ired = zeros(length(y),length(x));


True_Temp= zeros(length(x),length(y));

for i=1:maxnum
%This solves the non-dimensional temperature distribution for a rectangular plate when the
%top and bottom (T1) are at one temperature and the sides (T2) are at
%another temperature.  You can treat this as a black box.
for iX = 1:length(x)
    for iY = 1:length(y)
        sum = 0;
        for n=1:30
            term1 = ((-1)^(n+1) + 1) / n;
            term2 = sin(n*pi*x(iX) / L);
            term3 = sinh(n*pi*y(iY) / L);
            term4 = sinh(n*pi*W / L);
            sum = term1 * term2 * (term3/term4) + sum;
        end
        theta(iX,iY) = (2/pi) * sum;
%         %Uncomment this if you want to Add noise
%          theta(iX,iY) = theta(iX,iY) + randn/10;
    end
end


%%
%Makes it so that the temperature distribution is symmetrical
theta = (theta +fliplr(theta))/2;

%Converts nondimensional temperature to True Temperature
True_Temp(:,:) = theta * (T2-T1) + T1;


%%
fsize = 16;
figure;
 contourf(True_Temp);
 h = colorbar; set(get(h,'label'),'string','Temp (K)','FontWeight','bold','FontSize',fsize);
 set(gca,'FontWeight','bold','FontSize',fsize,'XTickLabel',{},'YTickLabel',{}); 
 saveas(gcf,'Temp Map - Temps.fig')
 
figure;
 contour(True_Temp);colormap(hot)
 h = colorbar; set(get(h,'label'),'string','Temp (K)','FontWeight','bold','FontSize',fsize);
 set(gca,'FontWeight','bold','FontSize',fsize,'XTickLabel',{},'YTickLabel',{});
 


%% 
% This creates a linear relationship between the temperature and pixel
% values of an image, and assumes a 16 bit pixel depth so that smaller
% temperature variations can be seen (a lot of greyscale cameras are 12
% bit, but you can change the bit depther easily
% look at Linear to Curve fit to NN comparison (Autosaved).xlsx


%Image 16 bit grayscale intensity due to temperature (290-380 range)
% m =     [-8.24	-54.03	-210.27	-272.26	-66.21 ;
%         3067.50	26378.00	105405.26	139946.00	36360.00];
%         
% %Image 16 bit grayscale intensity due to temperature (300-312 range)
m =     [-61.78	-405.22	-1577.00	-2041.90	-496.58 ;
        19213.00	132276.00	517539.07	673575.00	166133.00];

Izero   = uint16(zeros([size(True_Temp),3]));
Iblue   = uint16(m(1,1)*True_Temp + m(2,1));
Igreen  = uint16(m(1,2)*True_Temp + m(2,2));
Iyellow = uint16(m(1,3)*True_Temp + m(2,3));
Iorange = uint16(m(1,4)*True_Temp + m(2,4));
Ired    = uint16(m(1,5)*True_Temp + m(2,5));

%%
% close all
imBlue = Izero; imBlue(:,:,3) = Iblue*20;
imGreen = Izero; imGreen(:,:,2) = Igreen*2;
imYellow = Izero; imYellow(:,:,1) = 204*256;       imYellow(:,:,2) = 204*256;     imYellow(:,:,3) = Iyellow/1.5;
imOrange = Izero; imOrange(:,:,1) = 2^16;       imOrange(:,:,2) = Iorange/2;
imRed = Izero; imRed(:,:,1) = Ired*2;

figure;imshow(imBlue)
 saveas(gcf,'Temp Map - Blue.fig')
figure;imshow(imGreen)
 saveas(gcf,'Temp Map - Green.fig')
figure;imshow(imYellow)
 saveas(gcf,'Temp Map - Yellow.fig')
figure;imshow(imOrange)
 saveas(gcf,'Temp Map - Orange.fig')
figure;imshow(imRed)
 saveas(gcf,'Temp Map - Red.fig')

% imViolet = Izero; imOrange(:,:,1) = Iorange;    imOrange(:,:,3) = Iorange;
% figure;imshow(imViolet)
end
