%% apply NN to Exp data
clear
clc
FolderExp='C:\Users\colto\Desktop\research\TEMP_lab_home\Expanded theta experimentation';
FolderFig='C:\Users\colto\Desktop\research\TEMP_lab_home\Expanded theta experimentation\figure';
cd(FolderExp)
sampling=50;
freq_scan=freq_log(0.05,1,20);
freq_scan=freq_scan(2:end);
T_scan=305:10:305;
F_scan=2:-1:1;
theta0=[];
theta1=[];
for n1=1:length(T_scan)
    for n2=1:length(F_scan)
        name0=['T',num2str(T_scan(n1)),'_',num2str(F_scan(n2)),'_scan'];
        name_scan{n1,n2}=name0;
    end
end
WL800=load('WL800.txt');%% wavelengths of spectrometer
[Tn,Fn]=size(name_scan);
%%

for iTn=1:Tn;
    
    for iFn=1:Fn
        cd(FolderExp)
        close all
        clear Time FWHM Ratio theta
        clear SpectPeak ExpSpect Exp_data Peak_WL ExpSpectNorm
        clear ExpDataIP ExpDataNorm ExpSpect PumpInteg PumpPeak
        clear T_recon T_ACDC T_acdc
        clear exp_data Time PumpLaser ExpSpect PumpInteg PumpPeak SpectInteg
        clear SpectPeak1 SpectPI SpectIP ExpSpect_sort IS FWHM_s ExpSpectNorm
        clear Peak_WL SpectPeak3 T_out ExpData ExpDataNorm ExpDataNorm2 TimeDepen
        clear T_recon fft_T T_ACDC T_fit
        load Spectra_range
        exp_data=load([name_scan{iTn,iFn},'.txt']);
        cd(FolderFig)
        name0=name_scan{iTn,iFn};
        [e1,e2]=size(exp_data);
        s_input=Spectra_range+1;
        PumpLaser=exp_data(:,55:110);
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
        WL800=load('WL800.txt');%% wavelengths of spectrometer.
        WL=WL800(s_input);
        SpectPeak=zeros(size(exp_data,1),1);
        Peak_WL=zeros(size(exp_data,1),1);
        for isize=1:length(Time);
            [SpectPeak(isize,:),Peak_WL(isize,:)]=findpeak(WL,ExpSpect(isize,:),30);
            Ratio(isize,:)=SpectInteg(isize,:)/SpectPeak(isize,:);
            FWHM(isize,:)=fwhm(WL,ExpSpect(isize,:));
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
        %         FolderNN='\\MEETPC-0239\Data\Fluorescentie\fre_t_f_285_335\front_1mm\SPNN';
        %         FolderNN='\\MEETPC-0239\Data\Fluorescentie\fre_t_f_285_335\mid_3mm2\SPNN';
        FolderNN='C:\Users\colto\Desktop\research\TEMP_lab_home\Expanded theta experimentation\SPNN';
        cd(FolderNN) %% NN based on normalized spectrum and peak wavelength
        load input_Q.mat
        input_train=theta(:,input_Q);
        load param.nn
        shiftt=load('shift.nn');
        factorr=load('factor.nn');
        thetam_exp=apply_NN(param,shiftt,factorr,input_train);
        output_exp=thetam_exp(:,end);
        output_exp=output_exp';
        
        
        % %         FolderNN='\\MEETPC-0239\Data\Fluorescentie\fre_t_f_285_335\front_1mm\SPMat';
        %         FolderNN='\\MEETPC-0239\Data\Fluorescentie\fre_t_f_285_335\front_1mm\IPMat';
        %         clear net
        %         cd(FolderNN)
        %         load net.mat;
        %         load input_Q.mat
        %         input_train=theta(:,input_Q);
        %         output_exp=sim(net,input_train');
        cd(FolderFig)
        h7=figure(7);%% single NN plot
        [AX,H1,H2]=plotyy(Time,output_exp,Time,PumpInteg);grid on
        set(get(AX(1),'Ylabel'),'String','predicted temperature (K)','fontsize',15)
        set(get(AX(2),'Ylabel'),'String','pump intensity (counts)','fontsize',15)
        set(AX(2),'YLim',[0 max(PumpInteg)*1.05],'XLim',[0,max(Time)*1.001])
        set(AX(1),'YLim',[min(Ratio)-1 max(Ratio)+1],'XLim',[0,max(Time)*1.001])
        xlabel('Time (Sec)','Fontsize',15)
        set(AX(2),'YLim',[0 max(PumpInteg)*1.05],'XLim',[0,max(Time)*1.001])
        set(AX(1),'YLim',[min(output_exp)-0.1 max(output_exp)+0.1])
        saveas(h7,[name0,'_Temp'],'jpg')
        
        T_recon(:,1)=Time;T_recon(:,2)=PumpInteg;T_recon(:,3)=SpectInteg;
        T_recon(:,4)=SpectPeak;T_recon(:,5)=Peak_WL;T_recon(:,6)=FWHM;
        T_recon(:,7)=Ratio;T_recon(:,8)=output_exp;
        
        txt_temperature_name=[name0,'_TempData'];
        save(txt_temperature_name,'T_recon','-ascii','-tabs');
        
        p_weight=polyfit(Time,output_exp',3);T_dc=polyval(p_weight,Time);
        T_ac=output_exp'-T_dc;
        T_Ref=PumpInteg-mean(PumpInteg);
        
        % [params_acT,T_ac_fit] = sinefit(T_ac,tnew,freq_scan(iFn),0,0);
        T_ACDC(:,1)=Time;T_ACDC(:,2)=PumpInteg;T_ACDC(:,3)=T_dc;T_ACDC(:,4)=T_ac;
        txt_TACDC_name=[name0,'_ACDC'];
        save(txt_TACDC_name,'T_ACDC','-ascii','-tabs');
        
        h8=figure(8);
        subplot(211)
        plot(Time,output_exp,'.r',Time,T_dc,'-k','linewidth',2);
        ylim([min(output_exp)-2,max(output_exp)+2])
        xlabel('time (sec)','fontsize',12)
        ylabel('temperature (K)','fontsize',12)
        subplot(212)
        plot(Time,T_ac,'.r');
        ylim([min(T_ac)*1.1,max(T_ac)*1.1])
        xlabel('time (sec)','fontsize',12)
        ylabel('AC temperature (K)','fontsize',12)
        saveas(h8,[name0,'_ACDC'],'jpg')
        
        fsample=50;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [ampRef,phaRef,f]=fftNN(T_Ref,fsample);
        [amp,pha,f]=fftNN(T_ac,fsample);
        [Amp,indamp]=max(amp);
        Pha=pha(indamp);
        [AmpRef,indref]=max(ampRef);
        PhaRef=phaRef(indref);
        Freq_fft=f(indamp);
        if abs(Freq_fft-f(indref))>1
            Freq_fft=f(indref);
            Amp=amp(indref);
            Pha=pha(indref);
        end
        
        h9=figure(9);
        subplot(211)
        plot(f,ampRef,'-k','linewidth',2)
        xlim([f(indref-10),f(indref+100)])
        ylim([min(ampRef),max(ampRef)*1.05])
        xlabel('frequency (hz)','fontsize',12)
        ylabel('amplitude (counts)','fontsize',12)
        subplot(212)
        plot(f,amp,'-r','linewidth',2)
        xlim([f(indref-10),f(indref+100)])
        ylim([min(amp),max(amp)*1.1])
        xlabel('frequency (hz)','fontsize',12)
        ylabel('amplitude (K)','fontsize',12)
        saveas(h9,[name0,'_fft'],'jpg')
        
        fft_T(:,1)=f;fft_T(:,2)=amp;fft_T(:,3)=(pha);fft_T(:,4)=ampRef;fft_T(:,5)=(phaRef);
        txt_fft_name=[name0,'_FFT'];
        save(txt_fft_name,'fft_T','-ascii','-tabs')
        
        FAP_fft(iFn,1)=Freq_fft; FAP_fft(iFn,2)=Amp; FAP_fft(iFn,3)=Pha; FAP_fft(iFn,4)=PhaRef; FAP_fft(iFn,5)=Pha-PhaRef;
        FAP_name=['T_',num2str(T_scan(iTn)),'_FAP'];
        txt_FAP_namefft=[FAP_name,'_fft'];
        save(txt_FAP_namefft,'FAP_fft','-ascii','-tabs')
        Amp_fft(iTn,iFn)=Amp;Pha_fft(iTn,iFn)=Pha-PhaRef;PhaRef_fft(iTn,iFn)=PhaRef;
        
    end
    h10=figure(10);
    subplot(211);semilogx(FAP_fft(:,1),FAP_fft(:,2),'-*r');
    xlim([0.04,10])
    ylabel('Amplitude','Fontsize',12);xlabel('Frequency (Hz)','Fontsize',12)
    subplot(212);semilogx(FAP_fft(:,1),unwrap(FAP_fft(:,5))*180/pi,'-*r');
    xlim([0.04,10])
    ylabel('Phase (Deg)','Fontsize',10);xlabel('Frequency (Hz)','Fontsize',10)
    saveas(h10,txt_FAP_namefft,'jpg')
end
save Amp_fft.txt Amp_fft -ascii -tabs
save Pha_fft.txt Pha_fft -ascii -tabs
save PhaRef_fft.txt PhaRef_fft -ascii -tabs
h11=figure(11);
imagesc(freq_scan,T_scan,Amp_fft);
xlabel('Frequency (Hz)','Fontsize',12)
ylabel('Temperature (K)','Fontsize',12)
colorbar
saveas(h11,'Amp_all','jpg')
h12=figure(12);
imagesc(freq_scan,T_scan,unwrap(Pha_fft)*180/pi);
xlabel('Frequency (Hz)','Fontsize',12)
ylabel('Temperature (K)','Fontsize',12)
colorbar
saveas(h12,'Pha_fft','jpg')

%% capacity
% cd('\\MEETPC-0239\Data\Fluorescentie\fre_t_f_20130911\specNN_Matlab');
% cd('\\MEETPC-0239\Data\Fluorescentie\fre_t_f_20130914\cgNN2');
cd('\\MEETPC-0239\Data\Fluorescentie\fre_t_f_285_335\mid_3mm2\SCAN20130926\1MM\figure');
% cd('\\MEETPC-0239\Data\Fluorescentie\fre_t_f_20130911\IP_Matlab');
Amp_T=load('Amp_fft.txt');
Pha_T=load('Pha_fft.txt');
T_scan=305:10:305;
freq_scan=freq_log(0.05,1,20);
Amp_T_f=[freq_scan',Amp_T'];
Pha_T_f=[freq_scan',unwrap(Pha_T')*180/pi];
save Amp_T_f.txt Amp_T_f -ascii -tabs
save Pha_T_f.txt Pha_T_f -ascii -tabs

i=sqrt(-1);
T_omega=Amp_T.*exp(i*Pha_T);
omega=2*pi*freq_scan;
for iT=1:length(T_scan)
    C_omega(iT,:)=1./(T_omega(iT,:).*(i*omega));
    lestring{iT}=['T=',num2str(T_scan(iT)),'K'];
end
Amp_C=abs(C_omega);
Pha_C=unwrap(angle(C_omega))*180/pi;
real_C=real(C_omega);
imag_C=imag(C_omega);
save Real_cw.txt real_C -ascii -tabs
save Imag_cw.txt imag_C -ascii -tabs
% Pha_C=unwrap(angle(C_omega))*180/pi;
h13=figure(13);
semilogx(freq_scan,real(C_omega),'-.','linewidth',2)
legend(lestring)
xlim([0.045,1])
hold on
grid on
xlabel('Frequency (Hz)','Fontsize',12)
ylabel('real(c)','Fontsize',12)
saveas(h13,'Real_Cw','jpg')

h14=figure(14);
semilogx(freq_scan,imag_C,'-.','linewidth',2)
legend(lestring)
xlim([0.045,1])
hold on
grid on
xlabel('Frequency (Hz)','Fontsize',12)
ylabel('imag(c)','Fontsize',12)
saveas(h14,'Imag_Cw','jpg')

h15=figure(15);
loglog(freq_scan,abs(T_omega),'-.','linewidth',2)
legend(lestring)
xlim([0.045,1])
hold on
grid on
xlabel('Frequency (Hz)','Fontsize',12)
ylabel('amp','Fontsize',12)
saveas(h15,'abs_T','jpg')

h16=figure(16);
semilogx(freq_scan,unwrap(angle(T_omega))*180/pi,'-.','linewidth',2)
legend(lestring)
xlim([0.045,1])
grid on
hold on
xlabel('Frequency (Hz)','Fontsize',12)
ylabel('phase (degree)','Fontsize',12)
saveas(h16,'pha_T','jpg')

%% performing window function and fft aferward
cd('\\MEETPC-0239\Data\Fluorescentie\fre_t_f_285_335\mid_3mm2\SCAN20130926\3MM\figure');
sampling=50;
freq_scan=freq_log(0.05,1,20);
T_scan=305:10:305;
F_scan=20:-1:1;
theta0=[];
theta1=[];
for iTn=1:length(T_scan)
    for iFn=1:length(F_scan)
        clear fft_T
        name0=['T',num2str(T_scan(iTn)),'_',num2str(F_scan(iFn)),'_scan'];
        name_scan{iTn,iFn}=name0;
        data0=load([name0,'_ACDC']);
        t0=data0(:,1);
        pump=data0(:,2);
        pump_AC=pump-mean(pump);
        T_AC=data0(:,4);
        [ampT,phaT,f]=fftNN(T_AC,sampling);
        [Amp,indamp]=max(ampT);
        
        w1=hann(length(pump_AC));
        pump_AC_window=pump_AC.*w1;
        T_AC_window=T_AC.*w1;

        [ampRef,phaRef,f]=fftNN(pump_AC_window,sampling);
        [AmpRef,indref]=max(ampRef);
        PhaRef=phaRef(indref);
        [ampT2,phaT2,f]=fftNN(T_AC_window,sampling);
        [Amp2,indamp2]=max(ampT2);
        Pha=phaT2(indamp2);
        Freq_fft=f(indamp2);
        if abs(Freq_fft-f(indref))>1
            Freq_fft=f(indref);
            Amp=amp(indref);
            Pha=pha(indref);
        end
        
        h8=figure(8);
        subplot(211)
        plot(t0,pump_AC,'-k',t0,pump_AC_window,'-r','linewidth',2)
        xlabel('time (second)','fontsize',12)
        ylabel('pump (a.u.)','fontsize',12)
        subplot(212)
        plot(t0,T_AC,'-k',t0,T_AC_window,'-r','linewidth',2)
        xlabel('time (second)','fontsize',12)
        ylabel('T (K)','fontsize',12)
        saveas(h8,[name0,'_window'],'jpg')
        
        h9=figure(9);
        subplot(211)
        plot(f,ampRef,'-k','linewidth',2)
        xlim([f(indref-10),f(indref+100)])
        ylim([min(ampRef),max(ampRef)*1.05])
        xlabel('frequency (hz)','fontsize',12)
        ylabel('amplitude (counts)','fontsize',12)
        subplot(212)
        plot(f,ampT,'-r','linewidth',2)
        xlim([f(indref-10),f(indref+100)])
        ylim([min(ampT),max(ampT)*1.1])
        xlabel('frequency (hz)','fontsize',12)
        ylabel('amplitude (K)','fontsize',12)
        saveas(h9,[name0,'_fft'],'jpg')
        
        fft_T(:,1)=f;fft_T(:,2)=ampT;fft_T(:,3)=(phaT2);fft_T(:,4)=ampRef;fft_T(:,5)=(phaRef);
        txt_fft_name=[name0,'_FFT'];
        save(txt_fft_name,'fft_T','-ascii','-tabs')
        
        FAP_fft(iFn,1)=Freq_fft; FAP_fft(iFn,2)=Amp; FAP_fft(iFn,3)=Pha; FAP_fft(iFn,4)=PhaRef; FAP_fft(iFn,5)=Pha-PhaRef;
        FAP_name=['T_',num2str(T_scan(iTn)),'_FAP'];
        txt_FAP_namefft=[FAP_name,'_fft'];
        save(txt_FAP_namefft,'FAP_fft','-ascii','-tabs')
        Amp_fft(iTn,iFn)=Amp;Pha_fft(iTn,iFn)=Pha-PhaRef;PhaRef_fft(iTn,iFn)=PhaRef;
    end
end

%% %% extract full fft

% % cd('\\MEETPC-0239\Data\Fluorescentie\fre_t_f_20130911\specNN_Matlab');
% cd('\\MEETPC-0239\Data\Fluorescentie\fre_t_f_20130914\cgNN2');
% % cd('\\MEETPC-0239\Data\Fluorescentie\fre_t_f_20130911\IP_CG\');
% % cd('\\MEETPC-0239\Data\Fluorescentie\fre_t_f_20130911\IP_Matlab');
% freq_scan=freq_log(0.05,2,15);
% T_scan=240:10:320;
% F_scan=15:-1:1;
% for n1=3:3
%     for n2=1:length(F_scan)
%         name0=['T',num2str(T_scan(n1)),'_',num2str(F_scan(n2)),'_scan_FFT'];
%         name_scan{n1,n2}=name0;
%         data0=load(name0);
%         f=data0(:,1);
%         amp0=data0(:,2);
%         pha0=data0(:,3);
%         ampRef=data0(:,4);
%         phaRef=data0(:,5);
%         [ARef,index0]=hamonic_pickup(f,ampRef,freq_scan(n2));
%         A=amp0(index0);
%         P0=pha0(index0);
%         PRef=phaRef(index0);
%         P=P0-PRef;
%         T_hm=[f(index0),A,P,ARef',PRef];
%         save([name0,'_hm'],'T_hm','-ascii','-tabs')
%         saveas(h14,'Imag_Cw','jpg')
% h17=figure(17);
% loglog(f(index0),A,'-','linewidth',2)
% xlim([0.045,10])
% hold on
% grid on
% xlabel('Frequency (Hz)','Fontsize',12)
% ylabel('amp','Fontsize',12)
% saveas(h17,'Real_Cw','jpg')
%
% h18=figure(18);
% semilogx(f(index0),unwrap(P)*180/pi,'-','linewidth',2)
% xlim([0.045,10])
% grid on
% hold on
% xlabel('Frequency (Hz)','Fontsize',12)
% ylabel('phase','Fontsize',12)
%     end
% end