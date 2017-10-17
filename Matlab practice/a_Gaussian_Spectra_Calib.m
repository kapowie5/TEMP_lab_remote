%%  Used for simulation data
% V7 is meant for simulating spectra with typical uncertainties, and a
%       phase and amplitude decay based on a frequency scan
%V7_1 fixes the time spacing issue
% close all
clear all
clc
uncertainty_upload  %Current directory
% uncertainty_upload_laptop

%% Checking which source has 1% error
[val,index] = max(uVectorNums);
errorStr = uVectorNames{index};
if val < sum(uVectorNums)
    errorStr = 'uAll'
end

orig = pwd;

%% This section defines linear relationships as a function of temperature for the follow spectral features
%       -Peak Intensity (PI), Peak Wavelength (PWL), Full Width Half
%        Maximum (FWHM), and Baseline (0).
%These terms are given a letter to relate to the Gaussian distribution
%parameters
%       -Peak Intensity is A, Peak Wavelength is B, Full Width Half
%        Maximum is C, and Baseline is D.
%These linear relationships take the form A(T) = A_0 +S_A*T, where A is the
%       feature of interest as a function of temperature (T), A_0 is the
%       intercept, S_A is the slope or sensitivity of that feature to
%       temperature, and T is temperature.

%Spectral Features Intercepts - Based on Experimental Data
A_PI_0  = 70066;    %Intercept for Peak Intensity
B_PWL_0 = 530.57;      %Intercept for Peak Wavelength
C_FWHM_0= 12.51;        %Intercept for FWHM
D_0     = 1300;        %Intercept for Baseline

%Spectral Features Slopes - Loosely Based on Experimental Data
S_A     = -155.1;     %Slope for Peak Intensity
S_B     = 0.24459;          %Slope for Peak Wavelength
S_C     = 0.29668;          %Slope for FWHM
S_D     = 0;            %Slope for Baseline




%% This section defines various variables for the simulation
WL  = load('GreenSpectrometerWavelengths.txt');    %Vector of length nWL, of the wavelengths that will be simulated
WL_index= [1:length(WL)];
nWL     = length(WL);
T0_0    = 300;      %Intial temperature of the simulation before heating







%% Create Files
% [ Folder ] = SharedSpectraFolders;
Folder.Current = pwd;
% Folder.SpectraSaveParent=[Folder.Doc,Folder.mainDir,Folder.fspec];
Folder.SpectraSaveParent=Folder.Current;
str2 = '2017-10-05';
% str2 = input('Today (YYYY-MM-DD)','s');
str = [str2,', Freq Calibrated Simulated Spectras, 1% error at ',errorStr];
str3 = 'Oct 05';
cd(Folder.SpectraSaveParent)
mkdir(str);


%% Simulation section
D_0 = 20000;
for iT = 1:10
cd([Folder.SpectraSaveParent,'\',str])
dummy=[errorStr,'.txt']; dlmwrite(dummy,val);
    T0 = (iT-1)*2+T0_0;    
%     Tinf    = T0 + T_Rise;
               % iT
    
%     for iCalib = 1:10
        % Loop of calibration files
        clear SpecSave Spec T
            nt = 250;
            t = linspace(0,100,nt);
            
                %This loop will go through all the time steps.  This way we can
                %   modulate the temperature as a function of time
                for jt = 1:nt
                    
                    %Variable T is meant to basically superimpose the AC on the DC
                    %component
                    T_noise = ((randn/4)*0.01);
                    T(jt) =  T0 + T_noise;   
                    Ref(jt) = 0;

                    A_PI    = (A_PI_0) + (S_A)*T(jt);
                    B_PWL   = (B_PWL_0) + (S_B)*T(jt);
                    C_FWHM  = (C_FWHM_0) + (S_C)*T(jt);
                    D       = (D_0 ) + (S_D)*T(jt);     

                    %create gaussian spectra
                    for kWL = 1:nWL
                        Spec(jt,kWL) = (A_PI*exp(- ((WL(kWL)-B_PWL).^2) *(4*log(2)) / (C_FWHM)^2))+((randn/16)*uSN) + D + randn*4;
                    end       

                    curr = pwd;
                    %Create Green and IR Reference, then add to Spectra
                    % MAYBE REMOVE?
                    cd(orig);
                    GrnI = find(WL<535 & WL>529);
                    Spec_temp(1,:) = ((maketrap(WL_index,min(GrnI),min(GrnI)+2,max(GrnI)-2,max(GrnI))) * ((Ref(jt)+1)/2))*15000;% + uNoise;
                    Spec(jt,:) = Spec(jt,:) + Spec_temp(1,:);                       
                    cd(curr);
                    %Save new time vector, Temp, Resistance, and TC
                    tt(jt) = t(jt);
%                 end

                end
      
            calibName = ['T',num2str(T0),'K_calib.calib'];
            SpecSave = [tt' squeeze(Spec(:,:)) T'];
            dlmwrite(calibName,SpecSave);
            clear SpecSave   
            
            Spec_T(iT,:,:) = Spec;
iT
end
    T0
  
            


% %%
% figure;
% hold all
% for iiF=1:nWL
%     test(iiF)=Spec(Index_z0,1 );
% end
% plot(WL,test)
% for iiF=1:nWL
%     test(iiF)=Spec(Index_z0,2 );
% end
% plot(WL,test)
% for iiF=1:nWL
%     test(iiF)=Spec(Index_z0,5 );
% end
% plot(WL,test)
% for iiF=1:nWL
%     test(iiF)=Spec(Index_z0,10 );
% end
% plot(WL,test)
% for iiF=1:nWL
%     test(iiF)=Spec(Index_z0,50 );
% end
% plot(WL,test)
% for iiF=1:nWL
%     test(iiF)=Spec(Index_z0,100 );
% end
% plot(WL,test)
% xlabel('Wavelength (nm)','FontWeight','bold','FontSize', fsize)
% ylabel('Counts','FontWeight','bold','FontSize', fsize)
% set(gca,'FontWeight','bold','FontSize',11)
% saveas(gcf,'Fig - Select Spectras','png')
% saveas(gcf,'Fig - Select Spectras','fig')
% 
% %%
% figure;
% temp=squeeze(Spec(1,:,:));
% plot(WL,temp/1000)
% % clear temp
% xlabel('Wavelength (nm)','FontWeight','bold','FontSize', fsize)
% ylabel('Counts (thousands)','FontWeight','bold','FontSize', fsize)
% set(gca,'FontWeight','bold','FontSize', fsize)
% xlim([550,800])
% saveas(gcf,'Fig - All Spectras','png')
% saveas(gcf,'Fig - All Spectras','fig')
% 
% 
% %%
% figure;plot(t,Spec(1,:,600))
% xlabel('Time (s)','FontWeight','bold','FontSize', fsize)
% ylabel('Counts','FontWeight','bold','FontSize', fsize)
% title(['Intensity vs t at wavelength, ',num2str(WL(600)),'nm'])
% saveas(gcf,'Fig - I vs t','png')
% saveas(gcf,'Fig - I vs t','fig')
% 
% %%
% figure;plot(t,Spec(1,:,30))
% xlabel('Time (s)','FontWeight','bold','FontSize', fsize)
% ylabel('Counts','FontWeight','bold','FontSize', fsize)
% title('Intensity vs t for Green Reference')
% saveas(gcf,'Fig - Ref','png')
% saveas(gcf,'Fig - Ref','fig')
% 
% %%
% figure;plot(t,T(1,:))
% xlabel('Time (s)','FontWeight','bold','FontSize', fsize)
% ylabel('Temp (K)','FontWeight','bold','FontSize', fsize)
% title('Temperature')
% saveas(gcf,'Temperature','png')
% saveas(gcf,'Temperature','fig')
% 
% %%
% cd(Folder.Current)
% errorStr
% beep
% beep
% beep

% %%
% % cd('C:\Users\Troy\Dropbox\3 - Belgium Spider Silk\Programs\Troy''s spectra ones\')
% cd('C:\Users\User\Dropbox\3 - Belgium Spider Silk\Programs\Troy''s spectra ones\')
% LeuvenSpectralAnalysisProgram_V6freq_good
