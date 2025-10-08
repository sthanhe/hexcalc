% As matters stand at this point, this part must be programmed yourself. 
% In the following lines you will find an example programme 
% consisting of 3 modules which should make it easier to get started.

% Disclaimer: A lot of programming still needs to be done here, as there is a lot of hardcoding in this part.

clear

% The search path of the folder structure must be added to be able to call all functions
addpath("C:\Users\tbrunaue\Documents\MATLAB\sandTES_HEX\functions\")
addpath("C:\Users\tbrunaue\Documents\MATLAB\sandTES_HEX\")

number_Module = 3;
Module_name = ["SH" "EVAP" "ECO"];
filename = 'example_data_input.xlsx'; %File name of the exel file with the input data
storage_location = "\sandTES_HEX\EXAMPLE"; 
input_val = generate_input_tab(number_Module, Module_name, filename, storage_location); %generates a table from the Excel file with the input values

% Recuperator
% The objects must be created depending on the number of recuperators used for air preheating
rec1 = recuperator();
rec2 = recuperator();
Rec = [rec1 rec2];

%for the future assignment of the recuperators to the individual modules. 
% In this case, rec1 is assigned to the SH module and rec2 to the EVAP and ECO modules
rec_pos = boolean([1 0 0 ; 0 1 1]);

% Defines a constant value for the heat transfer on the outside of the pipe to the sand bed
al_a = 800; % W/m^2K

% Objects
% 3 objects are created for the 3 modules
SH = hex(al_a);
SH.name = "SH";
EVAP = hex(al_a);
EVAP.name = "EVAP";
ECO = hex(al_a);
ECO.name = "ECO";

% The created objects are assigned to a 1x3 object. 
% The same object is called with SH and Modules(1)
Modules = [SH EVAP ECO];

[Modules.flow_type] = deal("crosscurrent");
[Modules.process_medium] = deal("water");
Modules(1).rec_object = Rec;


EVAP.is_NU = true; % The evaporator is set to natural circulation
EVAP.h_Tr = EVAP.h_bed + 3; %The drum height is set at 3 metres above the bed 
EVAP.cn = 10; % The start circulation number is set to 10
P_Tr = 149.929 * 10^5; %Pa; The drum pressure is set


for n = 1:size(Modules,2)  
    Modules(n).da = input_val{"da",n+1};
    Modules(n).s = input_val{"s",n+1};
    Modules(n).K = input_val{"K",n+1};
    Modules(n).nCells = input_val{"nCells",n+1};
    Modules(n).mDot1 = input_val{"mDot1",n+1};
    Modules(n).mDot3 = input_val{"mDot3",n+1};
    Modules(n).nTube = input_val{"nTube",n+1};
    Modules(n).l = input_val{"l_Tube",n+1};
    Modules(n).L_bed = input_val{"l_bed",n+1};
    Modules(n).l_cavity = input_val{"l_cavity",n+1};
    Modules(n).w_channel = input_val{"w_channel",n+1};
    Modules(n).T1_start = input_val{"T1_start",n+1};
    Modules(n).T3_start = input_val{"T3_start",n+1};
    Modules(n).p_start = input_val{"p_start",n+1};
    Modules(n).Qdot_start = input_val{"Qdot_start",n+1};
end

% Startwerte
SH.TQ_calc(1);
SH.hex_calculation(1,0);

EVAP.mDot1 = SH.mDot1*EVAP.cn;
EVAP.T3_start = SH.T3right(end);
EVAP.TQ_calc(true);
EVAP.hex_calculation(1,0);

% HP_ECO.mDot1 = mDot1_neu;
ECO.T3_start = EVAP.T3right(end);
ECO.TQ_calc();
ECO.hex_calculation(1,0);


%Sandbed
[Modules.d_p] = deal(120 * 10^-6);
[Modules.h_bed] =  deal(3); %m 
[Modules.rho_p] =  deal(2650);
dp_bed = FluBed.deltaP(Modules(1).h_bed,Modules(1).eps_por,Modules(1).rho_p);
[Modules.p_bed] =  deal(10^5+dp_bed*0.7);

%compressor
comp = compressor();
comp.h_bed = Modules(1).h_bed;
comp.dp = Modules(1).d_p;
comp.rho_p = Modules(1).rho_p;
comp.calculate(1.7) %facor: dp = 1.7*dp_bed
[Modules.comp_object] = deal(comp);

% Rec
nRec = size(Rec,2);
for t = 1:nRec
    R = rec_pos(t,:);
    T3m = [Modules(R).T3];

    Th_in = mean(T3m);
    rec_d = recuperator();
    rec_d.Th_in = Th_in;
    rec_d.Tc_in = comp.Tout;

    Abed= sum(([Modules(R).L_bed] + [Modules(R).l_cavity]) .*[Modules(R).w_channel]);
    mf = hex.mflu(comp.dp,comp.rho_p,comp.p_bed,rec_d.Th_in,Abed,4);
    [Modules(R).mDot4] = deal(mf);
    rec_d.mDot4 = mf;
    rec_d.desing(comp.pout);

    Rec(t).recalculate(rec_d,mf,Th_in,comp.Tout,comp.pout);
    [Modules(R).T_flu] = deal(mf);

    len = sum([Modules(R).nCells]);
    pSand = NaN(1,len);
    pSand(:) = Modules(1).p_bed;
    ha1 = DryAir.h(comp.pout,Rec(t).Tc_out);
    ha2 = DryAir.h(pSand,T3m);
    Q4 = mf*(ha2-ha1)./len;

    start = 0;
    for i = find(R)
        Modules(i).QLoss_Rec = Q4(start+1:start+Modules(i).nCells);
        start = start + Modules(i).nCells;
    end
end

err=repmat(10,1,5);
count = 1;
X = NaN(1,50);
M = NaN(1,50);
while err(end)>1 || ~isConv
    hsatv=IF97.h(P_Tr,NaN,1);
    mDot1_neu = (abs(EVAP.QSum))/(hsatv-ECO.h1left(1));
    mDot1_alt = SH.mDot1;

    abw = abs((mDot1_neu-mDot1_alt)/mDot1_alt);
    check_abw = abw > 0.005;
    mDot1_neu(check_abw) = mDot1_alt+sign(mDot1_neu-mDot1_alt)*0.005*mDot1_alt;
  
    m_EVAP = mDot1_neu*EVAP.cn;

    % SH
    SH.mDot1 = mDot1_neu;
    SH.Qdot_start = abs(SH.QSum);
    SH.T1_start = IF97.T_sat(P_Tr);
    SH.p_start = P_Tr;
    SH.TQ_calc();
    SH.hex_calculation(1);

    %EVAP
    EVAP.mDot1 = m_EVAP;
    EVAP.Qdot_start = abs(EVAP.QSum);

    EVAP.T3_start = SH.T3right(end);
    dp_FR = hex.rho(P_Tr,EVAP.h1right(end))*9.81*EVAP.h_Tr;
    EVAP.p_start = P_Tr+dp_FR;
    EVAP.T1_start = IF97.T_sat(P_Tr)-0.1;
    EVAP.TQ_calc();
    EVAP.hex_calculation();
    pressure_diff_2 = EVAP.delta_p1(1)-dp_FR;

    % ECO
    ECO.mDot1 = mDot1_neu;
    ECO.Qdot_start = abs(ECO.QSum);
    ECO.T3_start = EVAP.T3right(end);
    ECO.TQ_calc();
    ECO.hex_calculation();

    P_Tr = ECO.p1(1);
    % Rec
    nRec = size(Rec,2);
    for t = 1:nRec
        R = rec_pos(t,:);
        T3m = [Modules(R).T3];

        Th_in = mean(T3m);
        rec_d = recuperator();
        rec_d.Th_in = Th_in;
        rec_d.Tc_in = comp.Tout;

        Abed= sum(([Modules(R).L_bed] + [Modules(R).l_cavity]) .*[Modules(R).w_channel]);
        mf = hex.mflu(comp.dp,comp.rho_p,comp.p_bed,rec_d.Th_in,Abed,4);
        [Modules(R).mDot4] = deal(mf);
        rec_d.mDot4 = mf;
        rec_d.desing(comp.pout);

        Rec(t).recalculate(rec_d,mf,Th_in,comp.Tout,comp.pout);
        [Modules(R).T_flu] = deal(mf);

        len = sum([Modules(R).nCells]);
        pSand = NaN(1,len);
        pSand(:) = Modules(1).p_bed;
        ha1 = DryAir.h(comp.pout,Rec(t).Tc_out);
        ha2 = DryAir.h(pSand,T3m);
        Q4 = mf*(ha2-ha1)./len;

        start = 0;
        for i = find(R)
            Modules(i).QLoss_Rec = Q4(start+1:start+Modules(i).nCells);
            start = start + Modules(i).nCells;
        end
    end

    error = [abs(SH.T1_start - EVAP.T1left(1))];
    M_error = mDot1_neu-mDot1_alt;

    M(count) = M_error;
    count = count + 1;

    err=circshift(err,-1);
    err(end)=mean(error);

    %convergence condition
    derr=diff(err);
    isConv=(all(derr<0) && err(end)/err(1)>0.7 && abs(M_error) < 0.5) || (max(abs(derr))<0.8 && abs(M_error) < 0.5);
    disp(err);
    disp(M_error)
end

output(Modules,"example_output",'C:\Users\tbrunaue\Documents\MATLAB\sandTES_HEX\EXAMPLE',2)
