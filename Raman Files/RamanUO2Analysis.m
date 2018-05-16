close all
% From "High temperature Raman study of UO2: A possible tool for in situ
% estimation of irradiation?induced heating" paper

filename = 'HT Raman UO2.txt';
delimiter = '\t';
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
fclose(fileID);
HTRamanUO2 = [dataArray{1:end-1}];
clearvars filename delimiter formatSpec fileID dataArray ans;
omega = HTRamanUO2(:,1);
HTRam = HTRamanUO2(:,2:end);
T = (20:10:650);%	30	40	50	60	70	80	90	100	110	120	130	140	150	160	170	180	190	200	210	220	230	240	250	260	270	280	290	300	310	320	330	340	350	360	370	380	390	400	410	420	430	440	450	460	470	480	490	500	510	520	530	540	550	560	570	580	590	600	610	620	630	640	650];

%%
close all
fsize = 16;
figure;
for i=1:64
    plot(omega,HTRam(:,i));
    hold all
end
xlabel('Wavenumber (cm^{-1})','FontWeight','bold','FontSize', fsize)
ylabel('Magnitude (a.u.)','FontWeight','bold','FontSize', fsize)
set(gca,'FontWeight','bold','FontSize', fsize)

figure;
plot(omega,HTRam(:,62)) %62
xlabel('Wavenumber (cm^{-1})','FontWeight','bold','FontSize', fsize)
ylabel('Magnitude (a.u.)','FontWeight','bold','FontSize', fsize)
set(gca,'FontWeight','bold','FontSize', fsize)

figure;
plot(T,max(HTRam(50:80,:)))
xlabel('Temp (^{o}C)','FontWeight','bold','FontSize', fsize)
ylabel('Max Magnitude (a.u.)','FontWeight','bold','FontSize', fsize)
set(gca,'FontWeight','bold','FontSize', fsize)

%%
figure
surf(omega,T,HTRam')
title('Original Sampling');

%% Interpolate 2D data field
X = repmat(omega,[1,64]);
Y = repmat(T,[575,1]);
[xq,yq] = meshgrid(325:1:1267,20:1:650);
vq = griddata(omega,T,HTRam',xq,yq);

figure;
mesh(xq,yq,vq)


%% create interpolated fit curve function that will take omega and T
[fitres,gof] = createFit(omega,T,HTRam);


