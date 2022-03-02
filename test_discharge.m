%dicharge HTX2
h_dis = hex();
h_dis.mDot1 = 1.39;
h_dis.mDot3 = 5.77;
h_dis.di = 12.8*10^-3;
h_dis.e = 3.6*10^-3;
h_dis.l = 66;
h_dis.N = 7;
h_dis.Ra = 0.5*10^-6;
h_dis.K = 0.04*10^-3;

Qdot1 = 1.31*10^6;
T1in = 273.15+59.97;
p1in = 23.5*10^6;
rho1in = IF97.v(p1in,T1in,NaN).^-1;
T3in = 273.15+328.86;

h_dis.TQ_diagram(Qdot1,T1in,rho1in,T3in);
h_dis.hex_calculation(2);
h_dis.Tx_diagram();
% h_dis.TQ_diagram(abs(h_dis.Qdot),h_dis.T1right(end),h_dis.rho1(end),h_dis.T3left(1))
% h_dis.px_diagram();