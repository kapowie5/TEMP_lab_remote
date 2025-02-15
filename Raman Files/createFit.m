function [fitresult, gof] = createFit(omega, T, HTRam)
%CREATEFIT(OMEGA,T,HTRAM)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : omega
%      Y Input : T
%      Z Output: HTRam
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 31-May-2017 16:15:59


%% Fit: 'untitled fit 1'.
[xData, yData, zData] = prepareSurfaceData( omega, T, HTRam );

% Set up fittype and options.
ft = 'linearinterp';

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft, 'Normalize', 'on' );

% Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, [xData, yData], zData );
% legend( h, 'untitled fit 1', 'HTRam vs. omega, T', 'Location', 'NorthEast' );
% % Label axes
% xlabel omega
% ylabel T
% zlabel HTRam
% grid on
% view( 38.1, 33.9 );


