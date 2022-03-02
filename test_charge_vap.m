%charge vap
% 
h_vap = hex();
h_vap.mDot1 = 1.25;
h_vap.mDot3 = 8;
h_vap.di = 12.8*10^-3;
h_vap.e = 3.6*10^-3;
h_vap.l = 66;
h_vap.N = 7;
h_vap.Ra = 0.5*10^-6;
h_vap.K = 0.04*10^-3;
Qdot1 = 3*10^6;
T1in = 273.15+600;
p1in = 15*10^6;
rho1in = IF97.v(p1in,T1in,NaN).^-1;
T3in = 273.15+50;

h_vap.TQ_diagram(Qdot1,T1in,rho1in,T3in);
h_vap.hex_calculation(1);
h_vap.Tx_diagram();
% h_vap.px_diagram();
% h_vap.TQ_diagram(abs(h_vap.Qdot),h_vap.T1left(1),h_vap.rho1(1),h_vap.T3right(end))