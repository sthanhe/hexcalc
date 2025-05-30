function [fitresult, gof] = fit_a(p, T, a)
%CREATEFIT(P,T,A)
%  Create a fit.
%
%  Data for 'a(p,T)' fit:
%      X Input: p
%      Y Input: T
%      Z Output: a
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 29-Mar-2023 13:07:07


%% Fit: 'a(p,T)'.
[xData, yData, zData] = prepareSurfaceData( p, T, a );

% Set up fittype and options.
ft = 'linearinterp';
opts = fitoptions( 'Method', 'LinearInterpolant' );
opts.ExtrapolationMethod = 'linear';
opts.Normalize = 'on';

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft, opts );

% Plot fit with data.
figure( 'Name', 'a(p,T)' );
h = plot( fitresult, [xData, yData], zData );
legend( h, 'a(p,T)', 'a vs. p, T', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'p', 'Interpreter', 'none' );
ylabel( 'T', 'Interpreter', 'none' );
zlabel( 'a', 'Interpreter', 'none' );
grid on
view( -28.7, 30.0 );


