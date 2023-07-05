clear
% HTX1
HTX1 = hex();
HTX1.name = "HTX1";
HTX1.flow_type = "crosscurrent";
HTX1.di = 
HTX1.e = 
HTX1.K = 
HTX1.nTube = 
HTX1.npass = 
HTX1.w_channel = 
HTX1.l = HTX1.w_channel*HTX1.npass;
HTX1.nCells = 100;
HTX1.process_medium = "water";
HTX1.mDot1 = 
HTX1.mDot3 = 
HTX1.L_bed = 

% %% HTX2
HTX2 = hex();
HTX2.name = "HTX2";
HTX2.flow_type = "crosscurrent";
HTX2.di =
HTX2.e = 
HTX2.K = 
HTX2.nTube = 
HTX2.npass = 
HTX2.w_channel = 
HTX2.l = HTX2.w_channel*HTX2.npass;
HTX2.nCells = 
HTX2.process_medium = "water";
HTX2.mDot1 = 
HTX2.mDot3 = 
HTX2.L_bed = 



%HTX
HTX = [HTX1 HTX2];
HTX(1).T1_start =
HTX(1).p_start= 
HTX(2).p_start= 
HTX(2).T3_start = 
HTX(1).Qdot_start = 
HTX(2).Qdot_start = 


HTX(1).h_Bett = 
HTX(1).d_p = 
HTX(1).p_bed = 
HTX(1).rho_p = 
HTX(1).eps_por = 
HTX(1).A_Bett = 
HTX(2).A_Bett =


HTX.TQ_diagram();
HTX.hex_calculation(1,1);
output(HTX,"Name",'Speicherort')


