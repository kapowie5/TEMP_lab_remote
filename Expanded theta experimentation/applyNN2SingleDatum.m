%% apply NN to Exp data
function predtemp = applyNN2SingleDatum(exp_data)
clear -regexp ?!exp_data
clc

% FolderExp='C:\Users\Troy\Dropbox\3 - Belgium Spider Silk\Programs\From Liwang\Neural Network thing\front\single test';
FolderExp=pwd;
FolderFig=FolderExp;
cd(FolderExp)
name_scan{1}='apply';
%%
%         close all

        load Spectra_range
       % for n1=1:16
        %    expdata0=load([name_scan{iTn,iFn},num2str(n1),'.txt']);
         %   exp_data = [exp_data;expdata0;];
        %end
        cd(FolderFig)
        s_input=Spectra_range+1;
        PumpLaser=exp_data(:,255:310);
        Time=exp_data(:,1)-exp_data(1,1);
        ExpSpect=exp_data(:,s_input+1);
        Spectra_base=min(ExpSpect,[],2)*ones(1,size(ExpSpect,2));
        ExpSpect=ExpSpect-Spectra_base;
        PumpInteg=sum(PumpLaser,2);
        
        SpectInteg=sum(ExpSpect,2);
        WL800=load('GreenSpectrometerWavelengths.txt');
        % WL800=WaveRead;%% wavelengths of spectrometer
        WL=WL800(s_input);
        SpectPeak=zeros(size(exp_data,1),1);
        Peak_WL=zeros(size(exp_data,1),1);
            [SpectPeak(1,:),Peak_WL(1,:)]=findpeak(WL,ExpSpect(1,:),30);
            Ratio(1,:)=SpectInteg(1,:)/SpectPeak(1,:);
            FWHM(1,:)=1;
        Norm_matrix=SpectPeak*ones(1,size(ExpSpect,2));
        ExpSpectNorm=ExpSpect./Norm_matrix;
        theta=[SpectPeak,SpectInteg,Ratio,FWHM,Peak_WL,ExpSpectNorm];
        
        
        
        
        
        
        %% choose proper NN
        %        FolderNNSP='\\MEETPC-0239\Data\Fluorescentie\newSan_1007\p1\SPNN';
%         FolderNNIP='\\MEETPC-0239\Data\Fluorescentie\newSan_1007\p1\IPNN';
        
        FolderNNIP=pwd;%,'\IPNN'];
        cd(FolderNNIP) %% NN based on IP
        load input_Q.mat
        input_trainIP=theta(:,input_Q);
        
        load param.nn
        shiftt=load('shift.nn');
        factorr=load('factor.nn');
        thetam_expIP=apply_NN(param,shiftt,factorr,input_trainIP);
        output_expIP=thetam_expIP(:,end);
        output_expIP=output_expIP';
        predtemp = output_expIP
        
        
        
    
