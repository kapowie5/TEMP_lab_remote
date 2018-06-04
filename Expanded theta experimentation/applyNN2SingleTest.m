%% apply NN to Exp data
function applyNN2SingleTest(varargin)
clear -regexp ?!varargin
clc

% FolderExp='C:\Users\Troy\Dropbox\3 - Belgium Spider Silk\Programs\From Liwang\Neural Network thing\front\single test';
FolderExp=pwd;
FolderFig=FolderExp;
cd(FolderExp)
sampling=50;
freq_scan=[0.2,0.5];
T_scan=300;
name_scan{1}='apply';
% name_scan{2}='T305_DCAC';
WL800=load('GreenSpectrometerWavelengths.txt');
% WL800=WaveRead;%% wavelengths of spectrometer
[Tn,Fn]=size(name_scan);
%%

for iTn=1:Tn;
    
    for iFn=1:Fn
        cd(FolderExp)
%         close all
        clear Time FWHM Ratio theta
        clear SpectPeak ExpSpect Exp_data Peak_WL ExpSpectNorm
        clear ExpDataIP ExpDataNorm ExpSpect PumpInteg PumpPeak
        clear T_recon T_ACDC T_acdc
        clear exp_data Time PumpLaser ExpSpect PumpInteg PumpPeak SpectInteg
        clear SpectPeak1 SpectPI SpectIP ExpSpect_sort IS FWHM_s ExpSpectNorm
        clear Peak_WL SpectPeak3 T_out ExpData ExpDataNorm ExpDataNorm2 TimeDepen
        clear T_recon fft_T T_ACDC T_fit
        load Spectra_range
        exp_data = [];
       % for n1=1:16
        %    expdata0=load([name_scan{iTn,iFn},num2str(n1),'.txt']);
         %   exp_data = [exp_data;expdata0;];
        %end
        if nargin==0
        exp_data = load('DD6.txt');
        applyNN2SingleDatum(exp_data(2,:))
        else
            exp_data = varargin{1};
        end
        
        h7=figure(7);grid on
        ylabel('predicted temperature (K)','fontsize',15)
        xlabel('Time (Sec)','Fontsize',15)
        legend('IP NN','SP NN')
 
        exp_data = exp_data(2:end,3:end-5);
        cd(FolderFig)
        name0=name_scan{iTn,iFn};
        [e1,e2]=size(exp_data);
        s_input=Spectra_range+1;
        PumpLaser=exp_data(:,255:310);
        Time=exp_data(:,1)-exp_data(1,1);
        ExpSpect=exp_data(:,s_input+1);
        Spectra_base=min(ExpSpect,[],2)*ones(1,size(ExpSpect,2));
        ExpSpect=ExpSpect-Spectra_base;
        PumpInteg=sum(PumpLaser,2);
        % [params_pump,PumpIntegFit] = sinefit(PumpInteg,Time,freq_scan(iFn),0,0);
        
        %         h1=figure(1);
        %         plot(Time,PumpInteg,'-b','linewidth',2)
        %         xlabel('time (Sec)','Fontsize',15)
        %         ylabel('pump intensity','Fontsize',15)
        %         set(gca,'fontsize',15)
        %         grid on
        %         saveas(h1,[name0,'_PumpIntens'],'jpg')
        SpectInteg=sum(ExpSpect,2);
        WL800=load('GreenSpectrometerWavelengths.txt');
        % WL800=WaveRead;%% wavelengths of spectrometer
        WL=WL800(s_input);
        SpectPeak=zeros(size(exp_data,1),1);
        Peak_WL=zeros(size(exp_data,1),1);
        for isize=1:length(Time);
            [SpectPeak(isize,:),Peak_WL(isize,:)]=findpeak(WL,ExpSpect(isize,:),30);
            Ratio(isize,:)=SpectInteg(isize,:)/SpectPeak(isize,:);
            FWHM(isize,:)=1;
        end
        Norm_matrix=SpectPeak*ones(1,size(ExpSpect,2));
        ExpSpectNorm=ExpSpect./Norm_matrix;
        UpOnes=ones(size(exp_data,1),1);
        T_out=UpOnes;
        theta=[SpectPeak,SpectInteg,Ratio,FWHM,Peak_WL,ExpSpectNorm];
        
        h2=figure(2);
        [AX,H1,H2]=plotyy(Time,SpectInteg,Time,PumpInteg);grid on
        set(get(AX(1),'Ylabel'),'String','integrated intensity (counts)','fontsize',15)
        set(get(AX(2),'Ylabel'),'String','pump intensity (counts)','fontsize',15)
        set(H1,'LineWidth',3);set(H2,'linewidth',3)
        set(AX(2),'YLim',[0 max(PumpInteg)*1.05],'XLim',[0,max(Time)*1.001])
        set(AX(1),'YLim',[min(SpectInteg)*0.99 max(SpectInteg)*1.005],'XLim',[0,max(Time)*1.001])
        set(AX(2),'YLim',[0 max(PumpInteg)*1.05])
        xlabel('time (sec)','Fontsize',15)
        saveas(h2,[name0,'_integI'],'jpg')
        
        h3=figure(3);
        [AX,H1,H2]=plotyy(Time,SpectPeak,Time,PumpInteg);grid on
        set(get(AX(1),'Ylabel'),'String','peak intensity (counts)','fontsize',15)
        set(get(AX(2),'Ylabel'),'String','pump intensity(counts)','fontsize',15)
        set(H1,'LineWidth',3);set(H2,'linewidth',3)
        set(AX(2),'YLim',[0 max(PumpInteg)*1.05],'XLim',[0,max(Time)*1.001])
        set(AX(1),'YLim',[min(SpectPeak)*0.99 max(SpectPeak)*1.005],'XLim',[0,max(Time)*1.001])
        xlabel('time (sec)','Fontsize',15)
        saveas(h3,[name0,'_peakI'],'jpg')
        
        h4=figure(4);
        [AX,H1,H2]=plotyy(Time,Peak_WL,Time,PumpInteg);grid on
        set(get(AX(2),'Ylabel'),'String','pump intensity (counts)','fontsize',15)
        set(get(AX(1),'Ylabel'),'String','peak wavelength (nm)','fontsize',15)
        set(AX(2),'YLim',[0 max(PumpInteg)*1.05],'XLim',[0,max(Time)*1.001])
        set(AX(1),'YLim',[min(Peak_WL)-0.1 max(Peak_WL)+0.1],'XLim',[0,max(Time)*1.001])
        xlabel('time (sec)','Fontsize',15)
        saveas(h4,[name0,'_PWL'],'jpg')
        
        h5=figure(5);
        [AX,H1,H2]=plotyy(Time,FWHM,Time,PumpInteg);grid on
        set(get(AX(2),'Ylabel'),'String','pump intensity (counts)','fontsize',15)
        set(get(AX(1),'Ylabel'),'String','fwhm (nm)','fontsize',15)
        set(AX(2),'YLim',[0 max(PumpInteg)*1.05],'XLim',[0,max(Time)*1.001])
        set(AX(1),'YLim',[min(FWHM)-1 max(FWHM)+1],'XLim',[0,max(Time)*1.001])
        xlabel('time (sec)','Fontsize',15)
        saveas(h5,[name0,'_FWHM'],'jpg')
        
        h6=figure(6);
        [AX,H1,H2]=plotyy(Time,Ratio,Time,PumpInteg);grid on
        set(get(AX(2),'Ylabel'),'String','pump intensity (counts)','fontsize',15)
        set(get(AX(1),'Ylabel'),'String','ratio (a.u.)','fontsize',15)
        set(AX(2),'YLim',[0 max(PumpInteg)*1.05],'XLim',[0,max(Time)*1.001])
        set(AX(1),'YLim',[min(Ratio)-1 max(Ratio)+1],'XLim',[0,max(Time)*1.001])
        xlabel('time (sec)','Fontsize',15)
        saveas(h6,[name0,'_Ratio'],'jpg')
        
        %% choose proper NN
        %        FolderNNSP='\\MEETPC-0239\Data\Fluorescentie\newSan_1007\p1\SPNN';
%         FolderNNIP='\\MEETPC-0239\Data\Fluorescentie\newSan_1007\p1\IPNN';
        
        FolderNNSP=pwd;
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
        cd(FolderNNSP) %% NN based on SP
        load input_Q.mat
        input_trainSP=theta(:,input_Q);
        load param.nn
        shiftt=load('shift.nn');
        factorr=load('factor.nn');
        thetam_expSP=apply_NN(param,shiftt,factorr,input_trainSP);
        output_expSP=thetam_expSP(:,end);
        output_expSP=output_expSP';
      
        cd(FolderFig)
        h7=figure(7);%% single NN plot
        plot(Time,output_expIP,'.k',Time,output_expSP,'.r');grid on
        ylabel('predicted temperature (K)','fontsize',15)
        xlabel('Time (Sec)','Fontsize',15)
        legend('IP NN','SP NN')
        saveas(h7,[name0,'_Temp'],'jpg')
        
        T_recon(:,1)=Time;T_recon(:,2)=PumpInteg;T_recon(:,3)=SpectInteg;
        T_recon(:,4)=SpectPeak;T_recon(:,5)=Peak_WL;T_recon(:,6)=FWHM;
        T_recon(:,7)=Ratio;T_recon(:,8)=output_expIP;T_recon(:,8)=output_expSP;
        sampleexp = sum(output_expIP(1000:1200))/201;
        txt_temperature_name=[name0,'_TempData'];
        save(txt_temperature_name,'T_recon','-ascii','-tabs');
        
        
        
    end
    
end
