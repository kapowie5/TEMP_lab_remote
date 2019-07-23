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
load('Temp Map - Temps (Circles)');

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
imGreen = Izero; imGreen(:,:,2) = Igreen;
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