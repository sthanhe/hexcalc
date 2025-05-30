function [fitresult, gof] = fit_s(p, T, s)
%CREATEFIT(P,T,S)
%  Create a fit.
%
%  Data for 's(p,T)' fit:
%      X Input: p
%      Y Input: T
%      Z Output: s
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 29-Mar-2023 13:05:36


%% Fit: 's(p,T)'.
[xData, yData, zData] = prepareSurfaceData( p, T, s );

% Set up fittype and options.
ft = 'linearinterp';
opts = fitoptions( 'Method', 'LinearInterpolant' );
opts.ExtrapolationMethod = 'linear';
opts.Normalize = 'on';

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft, opts );

% Plot fit with data.
figure( 'Name', 's(p,T)' );
h = plot( fitresult, [xData, yData], zData );
legend( h, 's(p,T)', 's vs. p, T', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'p', 'Interpreter', 'none' );
ylabel( 'T', 'Interpreter', 'none' );
zlabel( 's', 'Interpreter', 'none' );
grid on
view( 8.1, 28.8 );


