%charge con
h_con = hex();
h_con.mDot1 = 1.25;
h_con.mDot3 = 8;
h_con.di = 12.8*10^-3;
h_con.e = 3.6*10^-3;
h_con.l = 66;
h_con.N = 7;
h_con.Ra = 0.5*10^-6;
h_con.K = 0.04*10^-3;

Qdot1 = 3*10^6;
T1in = 273.15+600;
p1in = 15*10^6;
rho1in = IF97.v(p1in,T1in,NaN).^-1;
T3in = 273.15+50;

h_con.TQ_diagram(Qdot1,T1in,rho1in,T3in);
% h_con.hex_calculation(1);
% h_con.Tx_diagram();
% h_con.px_diagram();
% h_con.TQ_diagram(abs(h_con.Qdot),h_con.T1left(1),h_con.rho1(1),h_con.T3right(end))
