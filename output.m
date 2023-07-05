function output(object,name,location,type)
    if nargin < 4
        type = 1;
    end
    %Grafik
    object.plot_Tx(type);
    saveas(gcf,string(location)+'\Tx_'+name+'.jpg','jpg')
    object.plot_TQ(type);
    saveas(gcf,string(location)+'\TQ_'+name+'.jpg','jpg')

    % values
    hex_name = zeros(1,size(object,2));
    Twater_in = NaN(1,size(object,2));
    Twater_out = NaN(1,size(object,2));
    Tsand_in = NaN(1,size(object,2));
    Tsand_out = NaN(1,size(object,2));
    mDot_water = NaN(1,size(object,2));
    mDot_sand = NaN(1,size(object,2));
    mdot_water = NaN(1,size(object,2));
    mdot_sand = NaN(1,size(object,2));
    Qdot_water = NaN(1,size(object,2));
    Qdot_sand = NaN(1,size(object,2));
    Qdot_loss = NaN(1,size(object,2));
    alpha_in = NaN(1,size(object,2));
    alpha_out = NaN(1,size(object,2));
    k = NaN(1,size(object,2));
    p1in = NaN(1,size(object,2));
    p1out = NaN(1,size(object,2));
    delta_p = NaN(1,size(object,2));

    for n=1:size(object,2)
        hex_name(n) = object(n).name;
        if object(1).charge
            Twater_in(n) = object(n).T1left(1)-273.15;
            Twater_out(n) = object(n).T1right(end)-273.15;
            Tsand_in(n) = object(n).T3right(end)-273.15;
            Tsand_out(n) = object(n).T3left(1)-273.15;
        else
            Twater_in(n) = object(n).T1right(end)-273.15;
            Twater_out(n) = object(n).T1left(1)-273.15;
            Tsand_in(n) = object(n).T3left(1)-273.15;
            Tsand_out(n) = object(n).T3right(end)-273.15;
        end

        mDot_water(n) = object(n).mDot1(1);
        mDot_sand(n) = object(n).mDot3(1);
        mdot_water(n) = object(n).mDot1(1)/((object(n).di^2*pi()/4) * (object(n).nTube));
        mdot_sand(n) = object(n).mDot3(1)/(object(n).L_bed * object(n).w_channel);
        Qdot_water(n) = abs(sum(object(n).Qdot_1/10^6));
        Qdot_sand(n) = abs(sum(object(n).Qdot_3/10^6));
        Qdot_loss(n) = abs(sum(object(n).QLoss_Rec/10^6));
        alpha_in(n) =  mean(object(n).alpha_i);
        alpha_out(n) =  mean(object(n).alpha_a);
        k(n) = object(n).k;
        p1in(n) =  object(n).p1(end)/10^5;
        p1out(n) =  object(n).p1(1)/10^5;
        delta_p(n) = object(n).delta_p1(1)/10^5;
    end
    
    values_data = [hex_name; Twater_in; Twater_out; Tsand_in; Tsand_out; mDot_water; mDot_sand;...
        mdot_water; mdot_sand; Qdot_water; Qdot_sand; Qdot_loss; alpha_in; alpha_out;...
        k; p1in; p1out; delta_p];
    
    Name_data = ["name";"Twater_in"; "Twater_out"; "Tsand_in"; "Tsand_out"; "mDot_water";...
        "mDot_sand"; "mdot_water"; "mdot_sand"; "Qdot_water";
        "Qdot_sand"; "Qdot_loss"; "alpha_in"; "alpha_out"; "k"; "p1in"; "p1out";...
        "delta_p"];
        
    unit_data = ["-";"°C"; "°C"; "°C"; "°C"; "kg/s"; "kg/s"; "kg/sm^2"; "kg/sm^2"; "MW";  "MW"; "MW";...
        "W/m^2K"; "W/m^2K"; "W/m^2K"; "bar"; "bar"; "bar"];
   
    % compressor
    compressor = object(1).comp_object;
    Name_comp = ["dp"; "T_in"; "T_out"; "eta"; "wt"];
    unit_comp = ["bar"; "°C"; "°C"; "-"; "kJ/kg"];
    dp = zeros(1,size(object,2));
    T_in = zeros(1,size(object,2));
    T_out = zeros(1,size(object,2));
    eta = zeros(1,size(object,2));
    wt = zeros(1,size(object,2));

    dp(1) = compressor.dp_c/10^5;
    T_in(1) = compressor.Tin-273.15;
    T_out(1) = compressor.Tout-273.15;
    eta(1) = compressor.eta;
    wt(1) = compressor.wt/10^3;
  
    values_comp = [dp; T_in; T_out; eta; wt];

    % sand bed
    Name_sand = ["dp"; "rho"; "eps"; "dp_bed"; "fluid_degree"];
    unit_sand = ["µm"; "°kg/m3"; "-"; "bar"; "-"];
    dp = zeros(1,size(object,2));
    rho = zeros(1,size(object,2));
    eps = zeros(1,size(object,2));
    dp_bed = zeros(1,size(object,2));
    fluid_degree = zeros(1,size(object,2));
   
    dp(1) = object(1).d_p*10^6;
    rho(1) = object(1).rho_p;
    eps(1) = object(1).eps_por;
    dp_bed(1) = compressor.delta_p_bed/10^5;
    fluid_degree(1) = 4;
    
    values_sand= [dp; rho; eps; dp_bed; fluid_degree];


    % Recuperator
    rec = object(1).rec_object;
    Th_in2 = zeros(1,size(object,2));
    Tc_out2 = zeros(1,size(object,2));
    Tc_in2 = zeros(1,size(object,2));
    Th_out2 = zeros(1,size(object,2));
    QDot_rec2 = zeros(1,size(object,2));
    kA2 = zeros(1,size(object,2));
    m_fl2 = zeros(1,size(object,2));
    P_comp = zeros(1,size(object,2));
    for n=1:size(rec,2)
        Th_in2(n) = rec(n).Th_in-273.25;
        Tc_out2(n) = rec(n).Tc_out-273.25;
        Tc_in2(n) = rec(n).Tc_in-273.25;
        Th_out2(n) = rec(n).Th_out-273.25;
        QDot_rec2(n) = rec(n).QDot*20^-3;
        kA2(n) = rec(n).kA*20^-3;
        m_fl2(n) = rec(n).mDot4;
        P_comp(n) = m_fl2(n)*compressor.wt/10^3;
    end
    values_Rec = [Th_in2; Tc_out2; Tc_in2; Th_out2; QDot_rec2; kA2; m_fl2; P_comp];
    Name_Rec = ["TH_in"; "Tc_out"; "Tc_in"; "Th_out"; "QDot"; "kA"; "mDot_fl"; "P_comressor"];

    unit_Rec = ["°C"; "°C"; "°C"; "°C"; "kW"; "kW/K"; "kg/s"; "kW"];

    Name = [Name_data;Name_comp; Name_sand; Name_Rec];
    unit = [unit_data; unit_comp; unit_sand; unit_Rec];
    values = [values_data; values_comp; values_sand; values_Rec];
    T = table(Name,unit,values);
    file_name = string(location)+'\data_'+string(name)+".xlsx";
    writetable(T,file_name,'WriteRowNames',true)
end