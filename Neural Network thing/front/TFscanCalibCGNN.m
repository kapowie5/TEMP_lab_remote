%% calibrate NN for FandT scan.
clear
clc
close all
OHBOYNEWSTUFF!
FolderExp=pwd;
% FolderExp='\\MEETPC-0239\Data\Fluorescentie\newSan_1007\front';
cd(FolderExp)
sampling=250;
freq_scan=freq_log(0.05,2,15);
T_scan=295:10:335;
F_scan=15:-1:1;
theta0=[];
for n1=1:length(T_scan)
    for n2=1:length(F_scan)
        disp(T_scan(n1));
        name0=['T',num2str(T_scan(n1)),'_',num2str(F_scan(n2)),'_calib','.txt'];
        name_calib{n1,n2}=name0;
    end
    data0=load(['T',num2str(T_scan(n1)),'_calib','.txt']);
    theta0=[theta0;data0];
end
Spectra_range=[150:750];%[150:750];
save Spectra_range.mat Spectra_range
Spectra=theta0(:,3:end);
T_pt1000_all=theta0(:,1);
Spectra_all=theta0(:,3:end);
Spectra_base=min(Spectra_all,[],2)*ones(1,size(Spectra_all,2));
Spectra_all=Spectra_all-Spectra_base;
WL800=load('WL800.txt');%% wavelengths of spectrometer
[Tn,Fn]=size(name_calib);

%% calculate peak intensity, integrated intensity,fwhm, emission maximum
for ii=0:size(theta0,1)/sampling-1
    range_1=(1:sampling)+ii*sampling;
    theta_av(ii+1,:)=mean(theta0(range_1,:)); %
end
T_pt1000_av=theta_av(:,1);
Spectra_av=theta_av(:,3:end);
for isp=1:size(Spectra_av,1);
    [I_peak(isp,:),Peak_WL(isp,:)]=findpeak(WL800(Spectra_range),Spectra_av(isp,Spectra_range),30);
    I_integ(isp,:)=sum(Spectra_av(isp,Spectra_range));
    Ratio(isp,:)=I_integ(isp,:)/I_peak(isp,:);
    FWHM(isp,:)=fwhm(WL800(Spectra_range),Spectra_av(isp,Spectra_range));
    Spectra_norm_av(isp,:)=Spectra_av(isp,:)./I_peak(isp);
end

for isp=1:size(theta0,1)
    [I_peak_all(isp,:),Peak_WL_all(isp,:)]=findpeak(WL800(Spectra_range),Spectra_all(isp,Spectra_range),30);
    I_integ_all(isp,:)=sum(Spectra_all(isp,Spectra_range));
    Ratio_all(isp,:)=I_integ_all(isp,:)/I_peak_all(isp,:);
    FWHM_all(isp,:)=fwhm(WL800(Spectra_range),Spectra_all(isp,Spectra_range));
    Spectra_norm_all(isp,:)=Spectra_all(isp,:)./I_peak_all(isp);
end
theta_all=[T_pt1000_all,I_peak_all,I_integ_all,Ratio_all,FWHM_all,Peak_WL_all,Spectra_all(:,Spectra_range)];
theta_all_norm=[T_pt1000_all,I_peak_all,I_integ_all,Ratio_all,FWHM_all,Peak_WL_all,Spectra_norm_all(:,Spectra_range)];
save(['theta_all'],'theta_all','-ascii','-tabs');
save(['theta_all_norm'],'theta_all_norm','-ascii','-tabs');
% clear theta_all theta_all_norm Spectra_all theta0
%% plot temperature dependent spectra/normalized spectra
h11=figure(11);
box on
grid on
hold on

plot(WL800(Spectra_range),Spectra_av(:,Spectra_range),'linewidth',2);
xlabel('wavelength (nm)','fontsize',20)
ylabel('intensity (counts)','fontsize',20)
axis('tight') 
set(gca,'fontsize',15)
saveas(h11,'TempDepenSpecta','jpg')
h12=figure(12);
box on
grid on
hold on
plot(WL800(Spectra_range),Spectra_norm_av(:,Spectra_range),'linewidth',2);
xlabel('wavelength (nm)','fontsize',20)
ylabel('normalized intensity (a.u.)','fontsize',20)
axis('tight') 
set(gca,'fontsize',15)
saveas(h12,'TempDepenSpectaNorm','jpg')

%% temperature dependent intensity
h21=figure(21);
grid on
box on
hold on
plot(T_pt1000_av,I_integ,'-*b','linewidth',3);
xlabel('temperature (K)','fontsize',20)
ylabel('integrated intensity (counts)','fontsize',20)
set(gca,'fontsize',15)
saveas(h21,'TempDepenIntegI','jpg')

h22=figure(22);
grid on
box on
hold on
plot(T_pt1000_av,I_peak,'-*b','linewidth',3);
xlabel('temperature (K)','fontsize',20)
ylabel('peak intensity (counts)','fontsize',20)
set(gca,'fontsize',15)
saveas(h22,'TempDepenPeakI','jpg')

h23=figure(23);
box on
grid on
hold on
plot(T_pt1000_av,I_integ./I_peak,'-*b','linewidth',3);
xlabel('temperature (K)','fontsize',20)
ylabel('ratio (a.u.)','fontsize',20)
set(gca,'fontsize',15)
saveas(h23,'TempDepenRatio','jpg')

%% temperature dependent FWHM/ Emission maximum
h31=figure(31);
hold on
box on
grid on
plot(T_pt1000_av,FWHM,'-*b','linewidth',3);
xlabel('temperature (K)','fontsize',20)
ylabel('FWHM (nm)','fontsize',20)
set(gca,'fontsize',15)
saveas(h31,'TempDepenFWHM','jpg')

h32=figure(32);
box on
grid on
hold on
plot(T_pt1000_av,Peak_WL,'-*b','linewidth',3);
xlabel('temperature (K)','fontsize',20)
ylabel('emission maximum (nm)','fontsize',20)
set(gca,'fontsize',15)
saveas(h32,'TempDepenPeakWL','jpg')

%% save data for temperature dependent fluroescence
name='back_cali';
T_ipfm=[T_pt1000_av,I_peak,I_integ,Ratio,FWHM,Peak_WL];
save([name,'_TempDepen.txt'],'T_ipfm','-ascii','-tabs') 
TempDepenSpec=Spectra_av';
TempDepenSpecNorm=[T_pt1000_av,Spectra_norm_av]';
save([name,'_TempDepenSpec.txt'],'TempDepenSpec','-ascii','-tabs')
save([name,'_TempDepenSpecNorm.txt'],'TempDepenSpecNorm','-ascii','-tabs')

%% train NN
theta=load('theta_all_norm');
theta=[theta(:,2:end),theta(:,1)];%% [I,P,Ratio,FWHM,PWL,Norm]
theta=theta(1:end,:);%% [I,P,Ratio,FWHM,PWL,Norm]
%%
% input_Q=[1,2,5,round(linspace(6,606,60))];%%80
input_Q=[1:2];
save input_Q.mat input_Q
% load input_Q
factor4train=0.8;
SamplingF=250;
train_Q=[];test_Q=[];
% for itQ=1:size(theta,1)/SamplingF;
%     train_Q=[train_Q,(itQ-1)*SamplingF+1:itQ*SamplingF-round((1-factor4train)*SamplingF)];
%     test_Q=[test_Q,itQ*SamplingF-round((1-factor4train)*SamplingF)+1:SamplingF*itQ];
% end
% load train_Q.mat
% load test_Q.mat
RandIndex=randperm(size(theta,1));
train_Q=RandIndex(1:size(theta,1)*factor4train);
test_Q=RandIndex(size(theta,1)*factor4train+1:end);
% save train_Q.mat train_Q
% save test_Q.mat test_Q
input_train=theta(train_Q,input_Q);
target_train=theta(train_Q,end);
input_test=theta(test_Q,input_Q);
target_test=theta(test_Q,end);

%% Create Matlab Network
% numHiddenNeurons =4;  % Adjust as desired 10
% net = newfit(input_train',target_train',numHiddenNeurons);
% net.divideParam.trainRatio=90/100;  % Adjust as desired
% net.divideParam.valRatio=3/100;  % Adjust as desired
% net.divideParam.testRatio=2/100;  % Adjust as desired
% % Train and Apply Network
% [net,tr] = train(net,input_train',target_train');
% output_train=sim(net,input_train');
% output_test=sim(net,input_test');
% error_train=output_train'-target_train;
% error_test=output_test'-target_test;
% rms_train=(sum(error_train.*error_train)/length(error_train)).^(0.5)
% rms_test=(sum(error_test.*error_test)/length(error_test)).^(0.5)
% save net.mat net

%% train CG NN 
input_cg=[input_train,target_train;input_test,target_test];
target_T=input_cg(:,end);
cg1=size(input_train,1);
cg2=size(input_cg,2);
load param.nn
param([1,2])=[2,2];% configuration for NN
param([3,4])=[1,cg2-1];     % input colums
param([5,6])=[cg2,cg2];     % target colums
param([7,8])=[cg2,cg2];     % output colums
param(12)=50;               % iteration times
param([13,14])=[cg1,size(input_cg,1)]; % test rows,from param 13 to 14
param([15,16])=[2,1];       % hidden layers,output units
param(17)=1;       % norm
save param.nn param -ascii -tabs
thetam=train_NN2(param,input_cg);

output_T=thetam(:,end);
output_train=output_T(1:cg1);
output_test=output_T(cg1+1:end);
error_T=output_T-target_T;
error_train=output_train-target_train;
error_test=output_test-target_test;
rms_train=(sum(error_train.*error_train)/length(error_train)).^(0.5)
rms_test=(sum(error_test.*error_test)/length(error_test)).^(0.5)
load shift.nn
load factor.nn
thetam2=apply_NN(param,shift,factor,input_cg(:,1:end-1));
T_cg(:,1)=target_T;T_cg(:,2)=output_T;T_cg(:,3)=error_T;

%%
error_train=output_train-target_train;
error_test=output_test-target_test;
rms_train=(sum(error_train.*error_train)/length(error_train)).^(0.5)
rms_test=(sum(error_test.*error_test)/length(error_test)).^(0.5)
% save net_IPPW6_506_20.mat net
% Plot
% plotperf(tr)
%%
% close all
N_scan=floor(size(theta,1)/SamplingF);
[sortTrain,idTrain]=sort(target_train);
[sortTest,idTest]=sort(target_test);
T_av_train=ones(N_scan,2);
T_av_test=ones(N_scan,2);
for iN=0:N_scan-1
    range_1=(1:SamplingF)+iN*SamplingF;
    T_AV(iN+1,:)=[mean(theta(range_1,end))]; %
end
for iscan=1:length(T_AV)
    
idTrain0=find(target_train<(T_AV(iscan)+3)&(target_train>=(T_AV(iscan)-3)));
T_av_train(iscan,1)=mean(target_train(idTrain0));
T_av_train(iscan,2)=mean(output_train(idTrain0));

idTest0=find(target_test<(T_AV(iscan)+3)&(target_test>=(T_AV(iscan)-3)));
T_av_test(iscan,1)=mean(target_test(idTest0));
T_av_test(iscan,2)=mean(output_test(idTest0));

end
error_T_train_av=T_av_train(:,2)-T_av_train(:,1);
error_T_test_av=T_av_test(:,2)-T_av_train(:,1);

rms_train_av=(sum(error_T_train_av.*error_T_train_av)/length(error_T_train_av)).^(0.5);
rms_test_av=(sum(error_T_test_av.*error_T_test_av)/length(error_T_test_av)).^(0.5);

%%
% close all
figure(33)
subplot(211)
plot(target_train,error_train,'ob'); hold on
plot(target_test,error_test,'.r'); hold on
% xlim([260,325]);
ylabel('Absolute error (K)','Fontsize',10);xlabel('Target temperature (K)','Fontsize',10)
subplot(212)

plot(target_train,output_train,'ob'); hold on
plot(target_test,output_test,'.r'); hold on
name_rms_train_cg=['RMSE of train data: ',num2str(rms_train)];
name_rms_test_cg=['RMSE of test data: ',num2str(rms_test)];
text(290,320,name_rms_train_cg);text(290,310,name_rms_test_cg);
% ylim([215,325])
% xlim([215,325]);
ylabel('Output temperature (K)','Fontsize',10);xlabel('Target temperature (K)','Fontsize',10)
name33=['T_NN'];
savefig(name33)

figure(44)
subplot(211)
plot(T_av_train(:,2),error_T_train_av,'ob'); hold on
plot(T_av_train(:,1),error_T_test_av,'.r'); hold on
% xlim([215,325]);
ylabel('Absolute error (K)','Fontsize',10);xlabel('Target temperature (K)','Fontsize',10);
subplot(212)
plot(T_av_train(:,1),T_av_train(:,2),'ob'); hold on
plot(T_av_test(:,1),T_av_test(:,2),'.r'); hold on
name_rms_train_cg_av=['RMSE of train data: ',num2str(rms_train_av)];
name_rms_test_cg_av=['RMSE of test data: ',num2str(rms_test_av)];
text(290,320,name_rms_train_cg_av);text(290,310,name_rms_test_cg_av);
% ylim([215,325]);
% xlim([215,325]);
ylabel('Output temperature (K)','Fontsize',10);xlabel('Target temperature (K)','Fontsize',10)
name44=['T_NN_AV'];
savefig(name44)
error_T_train=[target_train,output_train,error_train];
error_T_test=[target_test,output_test,error_test];
error_T_trainAV=[T_av_train,error_T_train_av];
error_T_testAV=[T_av_test,error_T_test_av];
save error_T_train.txt error_T_train -ascii -tabs
save error_T_test.txt error_T_test -ascii -tabs
save error_T_trainAV.txt error_T_trainAV -ascii -tabs
save error_T_testAV.txt error_T_testAV -ascii -tabs

%%
applyNN2SingleTest
