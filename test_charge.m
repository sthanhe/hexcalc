%charge HTX2
h_ch = hex();
h_ch.mDot1 = 1.25;
h_ch.mDot3 = 4.81;
h_ch.di = 12.8*10^-3;
h_ch.e = 3.6*10^-3;
h_ch.l = 66;
h_ch.N = 7;
h_ch.Ra = 0.5*10^-6;
h_ch.K = 0.04*10^-3;

Qdot1 = 1.6*10^6;
T1in = 273.15+437.32;
p1in = 23.5*10^6;
rho1in = IF97.v(p1in,T1in,NaN).^-1;
T3in = 273.15+99.35;

h_ch.TQ_diagram(Qdot1,T1in,rho1in,T3in);
h_ch.hex_calculation(1);
h_ch.Tx_diagram();
% h_ch.px_diagram();
% h_ch.TQ_diagram(abs(h_ch.Qdot),h_ch.T1left(1),h_ch.rho1(1),h_ch.T3right(end))
