% % % % %dicharge HTX2
h_vapd = hex();
h_vapd.mDot1 = 1.39;
h_vapd.mDot3 = 9;
h_vapd.di = 12.8*10^-3;
h_vapd.e = 3.6*10^-3;
h_vapd.l = 66;
h_vapd.N = 8;
h_vapd.Ra = 0.5*10^-6;
h_vapd.K = 0.04*10^-3;
Qdot1 = 3*10^6;
T1in = 273.15+120;
p1in = 15*10^6;
rho1in = IF97.v(p1in,T1in,NaN).^-1;
T3in = 273.15+570;

h_vapd.TQ_diagram(Qdot1,T1in,rho1in,T3in);
h_vapd.hex_calculation(1);
h_vapd.Tx_diagram();
% h_vapd.px_diagram();
% h_vapd.TQ_diagram(abs(h_vapd.Qdot),h_vapd.T1right(end),h_vapd.rho1(end),h_vapd.T3left(1))