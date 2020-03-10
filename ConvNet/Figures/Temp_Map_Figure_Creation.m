close all
clear all

load('C:\Users\Munro\Desktop\singlepixelgeneralizationfigures (1)\singlepixelgeneralizationfigures\30x30ErrorMatrix.mat')
load('C:\Users\Munro\Desktop\singlepixelgeneralizationfigures (1)\singlepixelgeneralizationfigures\30x30PredMatrix.mat')
load('C:\Users\Munro\Desktop\singlepixelgeneralizationfigures (1)\singlepixelgeneralizationfigures\30x30TrueMatrix.mat')

%%
figure1 = figure;

%Input Figure
subplot(3,3,1)
contour(inputmat)
set(gca,'BoxStyle','full','Layer','top','XTick',zeros(1,0),'YTick',...
     zeros(1,0));
title('Input','FontSize',14,'FontWeight','bold')
ylabel({'Pixel(x,y)';'(30,30)'},'FontSize',12,'FontWeight','bold')
colorbar(subplot1,'Ticks',[380 385 390],'TickLabels',{'380','385','390'},...
    'Limits',[380 390],...
    'FontWeight','bold',...
    'FontSize',12);

%Output Figure
% figure;
subplot(3,3,2)
contour(outputmat)
set(gca,'BoxStyle','full','Layer','top','XTick',zeros(1,0),'YTick',...
     zeros(1,0));
title('NN Output','FontSize',14,'FontWeight','bold')
% axes1 = axes(gca)

%Error Figure
% figure;
subplot(3,3,3)
contour(Errormat);
set(gca,'BoxStyle','full','Layer','top','XTick',zeros(1,0),'YTick',...
     zeros(1,0));
title('Error','FontSize',14,'FontWeight','bold')
 
 
 

load('C:\Users\Munro\Desktop\singlepixelgeneralizationfigures (1)\singlepixelgeneralizationfigures\60x100ErrorMatrix.mat')
load('C:\Users\Munro\Desktop\singlepixelgeneralizationfigures (1)\singlepixelgeneralizationfigures\60x100PredMatrix.mat')
load('C:\Users\Munro\Desktop\singlepixelgeneralizationfigures (1)\singlepixelgeneralizationfigures\60x100TrueMatrix.mat')



%Input Figure
subplot(3,3,4)
contour(inputmat)
set(gca,'BoxStyle','full','Layer','top','XTick',zeros(1,0),'YTick',...
     zeros(1,0));
ylabel({'Pixel(x,y)';'(60,100)'},'FontSize',12,'FontWeight','bold')

%Output Figure
% figure;
subplot(3,3,5)
contour(outputmat)
set(gca,'BoxStyle','full','Layer','top','XTick',zeros(1,0),'YTick',...
     zeros(1,0));
% axes1 = axes(gca)

%Error Figure
% figure;
subplot(3,3,6)
contour(Errormat);
set(gca,'BoxStyle','full','Layer','top','XTick',zeros(1,0),'YTick',...
     zeros(1,0));




 
%%
load('C:\Users\Munro\Desktop\singlepixelgeneralizationfigures (1)\singlepixelgeneralizationfigures\190x110ErrorMatrix.mat')
load('C:\Users\Munro\Desktop\singlepixelgeneralizationfigures (1)\singlepixelgeneralizationfigures\190x110PredMatrix.mat')
load('C:\Users\Munro\Desktop\singlepixelgeneralizationfigures (1)\singlepixelgeneralizationfigures\190x110TrueMatrix.mat')

%%

%Input Figure
subplot(3,3,7)
contour(inputmat)
set(gca,'BoxStyle','full','Layer','top','XTick',zeros(1,0),'YTick',...
     zeros(1,0));
ylabel({'Pixel(x,y)';'(190,110)'},'FontSize',12,'FontWeight','bold')

%Output Figure
% figure;
subplot(3,3,8)
contour(outputmat)
set(gca,'BoxStyle','full','Layer','top','XTick',zeros(1,0),'YTick',...
     zeros(1,0));
% axes1 = axes(gca)

%Error Figure
% figure;
subplot(3,3,9)
contour(Errormat);
set(gca,'BoxStyle','full','Layer','top','XTick',zeros(1,0),'YTick',...
     zeros(1,0));

 
 
 
 
%%
load('C:\Users\Munro\Desktop\singlepixelgeneralizationfigures (1)\singlepixelgeneralizationfigures\circleErrorMatrix.mat')
load('C:\Users\Munro\Desktop\singlepixelgeneralizationfigures (1)\singlepixelgeneralizationfigures\circlePredMatrix.mat')
load('C:\Users\Munro\Desktop\singlepixelgeneralizationfigures (1)\singlepixelgeneralizationfigures\circleTrueMatrix.mat')

figure;
% %Input Figure
% subplot(1,2,1)
% contour(inputmat)
% set(gca,'BoxStyle','full','Layer','top','XTick',zeros(1,0),'YTick',...
%      zeros(1,0));
% title('Input','FontSize',14,'FontWeight','bold')
% ylabel({'Pixel(x,y)';'(30,30)'},'FontSize',12,'FontWeight','bold')

%Output Figure
subplot(2,2,1)
contour(outputmat)
set(gca,'BoxStyle','full','Layer','top','XTick',zeros(1,0),'YTick',...
     zeros(1,0));
title('NN Output','FontSize',14,'FontWeight','bold')
% axes1 = axes(gca)

%Error Figure
subplot(2,2,2)
contour(Errormat);
set(gca,'BoxStyle','full','Layer','top','XTick',zeros(1,0),'YTick',...
     zeros(1,0));
title('Error','FontSize',14,'FontWeight','bold')






 

load('C:\Users\Munro\Desktop\singlepixelgeneralizationfigures (1)\singlepixelgeneralizationfigures\rvszErrorMatrix.mat')
load('C:\Users\Munro\Desktop\singlepixelgeneralizationfigures (1)\singlepixelgeneralizationfigures\rvszPredMatrix.mat')
load('C:\Users\Munro\Desktop\singlepixelgeneralizationfigures (1)\singlepixelgeneralizationfigures\rvszTrueMatrix.mat')

% figure;
% %Input Figure
% subplot(1,2,1)
% contour(inputmat)
% set(gca,'BoxStyle','full','Layer','top','XTick',zeros(1,0),'YTick',...
%      zeros(1,0));
% title('Input','FontSize',14,'FontWeight','bold')
% ylabel({'Pixel(x,y)';'(30,30)'},'FontSize',12,'FontWeight','bold')

%Output Figure
subplot(2,2,3)
contour(outputmat)
set(gca,'BoxStyle','full','Layer','top','XTick',zeros(1,0),'YTick',...
     zeros(1,0));
title('NN Output','FontSize',14,'FontWeight','bold')
% axes1 = axes(gca)

%Error Figure
subplot(2,2,4)
contour(Errormat);
set(gca,'BoxStyle','full','Layer','top','XTick',zeros(1,0),'YTick',...
     zeros(1,0));
title('Error','FontSize',14,'FontWeight','bold')




 
% 
% % Create axes
% axes1 = axes('Parent',figure1,...
%     'Position',[0.0764285714285713 0.709264705882353 0.213405797101449 0.215735294117647]);
% hold(axes1,'on');
% box(axes1,'on');
% axis(axes1,'tight');
% set(axes1,'BoxStyle','full','Layer','top','XTick',zeros(1,0),'YTick',...
%     zeros(1,0));
% 
% % Create axes
% axes2 = axes('Parent',figure1,...
%     'Position',[0.319725672877846 0.709264705882353 0.210631469979296 0.215735294117647]);
% hold(axes2,'on');
% 
% % 
% box(axes2,'on');
% axis(axes2,'tight');
% % Set the remaining axes properties
% set(axes2,'BoxStyle','full','Layer','top','XTick',zeros(1,0),'YTick',...
%     zeros(1,0));
% % Create colorbar
% colorbar(axes2);
% 
% % Create axes
% axes3 = axes('Parent',figure1,...
%     'Position',[0.688022774327119 0.711645658263306 0.203048654244308 0.215735294117647]);
% hold(axes3,'on');
% 
% box(axes3,'on');
% axis(axes3,'tight');
% % Set the remaining axes properties
% set(axes3,'BoxStyle','full','Layer','top','XTick',zeros(1,0),'YTick',...
%     zeros(1,0));
% % Create colorbar
% colorbar(axes3);