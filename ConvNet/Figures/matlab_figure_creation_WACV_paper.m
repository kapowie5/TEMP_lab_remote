% matlab_figure_creation_WACV_paper
close all
fsize = 16;
msize = 8;
lsize = 1;

load('Spectral_curves_high_and_low')

% Spectral fits
figure;
semilogy(T_high,specB,'bo','MarkerFaceColor','b','MarkerSize',msize)
hold all
semilogy(T_high,specG,'gV','MarkerFaceColor','g','MarkerSize',msize)
semilogy(T_high,specY,'y^','MarkerFaceColor','y','MarkerSize',msize)
semilogy(T_high,specO,'d','MarkerEdgeColor',[0.9,0.41,0.17],'MarkerFaceColor',[0.9,0.41,0.17],'MarkerSize',msize)
semilogy(T_high,specR,'rs','MarkerFaceColor','r','MarkerSize',msize)

semilogy(T_high,fitB_spec,'b--','LineWidth',lsize)
semilogy(T_high,fitG_spec,'g--','LineWidth',lsize)
semilogy(T_high,fitY_spec,'y--','LineWidth',lsize)
semilogy(T_high,fitO_spec,'--','Color',[0.9,0.41,0.17],'LineWidth',lsize)
semilogy(T_high,fitR_spec,'r--','LineWidth',lsize)


set(gca,'YMinorTick','on','FontSize',fsize,'FontWeight','bold')
ylabel('Summed Intensity (counts)','FontSize',fsize,'FontWeight','bold')
xlabel('Temp (K)','FontSize',fsize,'FontWeight','bold')
xlim([280,390])
%%
%Grayscale fits
figure;
semilogy(T_high,grayB,'bo','MarkerFaceColor','b','MarkerSize',msize)
hold all
semilogy(T_high,grayG,'gV','MarkerFaceColor','g','MarkerSize',msize)
semilogy(T_high,grayY,'y^','MarkerFaceColor','y','MarkerSize',msize)
semilogy(T_high,grayO,'d','MarkerEdgeColor',[0.9,0.41,0.17],'MarkerFaceColor',[0.9,0.41,0.17],'MarkerSize',msize)
semilogy(T_high,grayR,'rs','MarkerFaceColor','r','MarkerSize',msize)

semilogy(T_high,fitB_gray,'b--','LineWidth',lsize)
semilogy(T_high,fitG_gray,'g--','LineWidth',lsize)
semilogy(T_high,fitY_gray,'y--','LineWidth',lsize)
semilogy(T_high,fitO_gray,'--','Color',[0.9,0.41,0.17],'LineWidth',lsize)
semilogy(T_high,fitR_gray,'r--','LineWidth',lsize)


set(gca,'YMinorTick','on','FontSize',fsize,'FontWeight','bold')
ylabel('Grayscale Intensity','FontSize',fsize,'FontWeight','bold')
xlabel('Temp (K)','FontSize',fsize,'FontWeight','bold')
xlim([280,390])