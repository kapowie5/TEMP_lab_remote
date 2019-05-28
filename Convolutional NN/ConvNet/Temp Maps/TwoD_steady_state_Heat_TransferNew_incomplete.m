clear all
close all

%Vary These temperatures to make different images
T1 = 20+273-5; %Top and bottom temperatures (K)
T2 = 96+273-5; %Side Temperatures (K)
%Min and max temps,   and  20+273

%Dimensions of your image
L = 121;  %pixels
W = 201;  %pixels


x = 0:1:L-1;
y = 0:1:W-1;
maxnum = 1;
Igreen = zeros(length(y),length(x));

Iblue = zeros(length(y),length(x));

Ired = zeros(length(y),length(x));


True_Temp= zeros(maxnum,length(x),length(y));

for i=1:maxnum
T1=T1+1;
T2=T2+2;
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
         theta(iX,iY) = theta(iX,iY) + randn/10;
    end
end


%%
%Makes it so that the temperature distribution is symmetrical
theta = (theta +fliplr(theta))/2;

%Converts nondimensional temperature to True Temperature
True_Temp(i,:,:) = theta * (T2-T1) + T1;
%% to create a circle (page 103 of Heat conduction, Latif M. Jiji, solution 108)


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
 True_Temp = T_full;
%%
% d_theta = 360/nr;
% theta  = [1:d_theta:360];

% [x,y] = pol2cart(theta,r);

%% 
% This creates a linear relationship between the temperature and pixel
% values of an image, and assumes a 16 bit pixel depth so that smaller
% temperature variations can be seen (a lot of greyscale cameras are 12
% bit, but you can change the bit depther easily
bit = 16;   %number of bits per pixel
grey_max = 2^bit - 1;
grey_min = 0;
temp_max = T2;
temp_min = T1;

m = (grey_max-grey_min) / (temp_max-temp_min);
b = grey_min - (m*temp_min);

%Image intensity due to temperature
I = m * True_Temp + b;
Ired(:,:) = -287.74*True_Temp(:,:) + 144236;
Igreen(:,:) = -73.933*True_Temp(:,:) + 36096;
Iblue(:,:) = -11.272*True_Temp(:,:) + 4198;
%%
end
test1 = squeeze(Ired(1,:,:));
test2 = squeeze(Ired(2,:,:));

TempList = zeros(201*121,1);
RedList = zeros(1,201*121,1);
GreenList = zeros(1,201*121,1);
BlueList = zeros(1,201*121,1);

for xl=1:201
   for yl = 1:121
       TempList(xl*yl)=True_Temp(xl,yl);
       RedList(xl*yl)=Ired(xl,yl);
       GreenList(xl*yl)=Igreen(xl,yl);
       BlueList(xl*yl)=Iblue(xl,yl);
   end
end
%Converts the intensity to a 16 bit integer for saving as an image
%I2 = uint16(I);


%imwrite(I2, '16bit3.png', 'BitDepth',16);

%figure;
%contour(y,x,True_Temp)

%figure;
%contour(y,x,I)

%figure;
%contour(y,x,I2)

%%
%Converts grayscale to RGB
%I3 = im2uint8(I2);

%I4 = grs2rgb(I3); %https://www.mathworks.com/matlabcentral/fileexchange/13312-grayscale-to-rgb-converter
%figure;imshow(I4)
%axes1 = axes('Parent',gcf,...
 %   'Position',[0.0586776845095572 0.0801635972073717 0.835629017447199 0.883435582822086]);
%axis off
%hold(axes1,'on');
%box(axes1,'on');
%axis(axes1,'ij');
% Set the remaining axes properties
%set(axes1,'CLim',[290 370],'DataAspectRatio',[1 1 1],'FontSize',16,...
 %   'FontWeight','bold','Layer','top','TickDir','out');
% Create colorbar
%colorbar('peer',axes1);

%set(gca,'FontSize',16,'FontWeight','bold')
%h = colorbar;
%ylabel(h, 'Temp (K)','FontSize',16,'FontWeight','bold')



%%
%Displays red, green, and blue channels
%I4 = grs2rgb(I3,hsv); %https://www.mathworks.com/matlabcentral/fileexchange/13312-grayscale-to-rgb-converter
%figure;imshow(I4)

%img = I4;
%red = img(:,:,1); % Red channel
%green = img(:,:,2); % Green channel
%blue = img(:,:,3); % Blue channel
%a = zeros(size(img, 1), size(img, 2));
%just_red = cat(3, red, a, a);
%just_green = cat(3, a, green, a);
%just_blue = cat(3, a, a, blue);
%back_to_original_img = cat(3, red, green, blue);
%figure, imshow(img), title('Original image')
%figure, imshow(just_red), title('Red channel')
%figure, imshow(just_green), title('Green channel')
%figure, imshow(just_blue), title('Blue channel')
%figure, imshow(back_to_original_img), title('Back to original image')
%figure;imshow(I4(:,:,1))
%figure;imshow(I4(:,:,2))
%figure;imshow(I4(:,:,3))

%% Write file with matrix of what temperatures should be
%csvwrite(['True_Temp_T1=',num2str(T1),'_T2=',num2str(T2)],True_Temp);


%% Random noise with some temp highlights
%I5 = zeros(1080,1920,3);
%I5(:,:,1) = imnoise(I5(:,:,1),'gaussian',0.1,0.05);
%I5(:,:,2) = imnoise(I5(:,:,2),'gaussian',0.1,0.002);
%I5(:,:,3) = imnoise(I5(:,:,3),'gaussian',0.1,0.01);
%Red square
%r=20;x=3;y=6;
%th = 0:pi/50:2*pi;
%xunit = r * cos(th) + x;
%yunit = r * sin(th) + y;
% I5(:,:,1) = imnoise(I5(xunit+500,yunit+1000,1),'gaussian',0.1,0.05);
%I5(xunit+500,yunit+1000,3) = I5(xunit+500,yunit+1000,3) + rand/2+0.5;
%I5(xunit+710,yunit+100,3) = I5(xunit+710,yunit+100,3) + rand/2+0.5;
%I5(xunit+ceil(rand*400),yunit+ceil(rand*400),3) = I5(xunit+ceil(rand*400),yunit+ceil(rand*400),3) + rand/2+0.5;
% I5 = insertShape(I5,'circle',[150 280 35],'LineWidth',5);
%figure;imshow(I5)

%%
%img = I5;
%red = img(:,:,2); % Red channel
%green = img(:,:,3); % Green channel
%blue = img(:,:,1); % Blue channel
%a = zeros(size(img, 1), size(img, 2));
%just_red = cat(3, red, a, a);
%just_green = cat(3, a, green, a);
%just_blue = cat(3, a, a, blue);
%back_to_original_img = cat(3, red, green, blue);
%figure, imshow(img), title('Original image')
%figure, imshow(just_red), title('Red channel')
%figure, imshow(just_green), title('Green channel')
%figure, imshow(just_blue), title('Blue channel')
%figure, imshow(back_to_original_img), title('Back to original image')
%figure;imshow(I4(:,:,1))
%figure;imshow(I4(:,:,2))
%figure;imshow(I4(:,:,3))

%%
%I6 = rgb2gray(I5);
%I6(I6<0.20) = 1;
%I7 = grs2rgb(I6);
%figure;imshow(I7)