classdef hex < handle 
    %RECALCULATION HEAT EXCAHNGER
    %references: VDI Wärmeatlas, 11.Auflage, SpringerVieweg
                  
     % charge/discharge
     %1: H2O-Side
     %3: Sand-Side
     
    properties
        name
        flow_type
        process_medium
        prop = IF97;
        % start conditions
        Qdot_start = 0
        T1_start = 0
        p_start = 0
        T3_start = 0

        mDot1;       % mass flow process medium(kg/s)                 [1,1]
        mDot3;       % mass flow sand (kg/s)                          [1,1]
        mDot4;       % mass flow fluidisation air (kg/s)
        di;          % inner diameter (m)                             [1,1]
        e;           % wall thickness (m)                             [1,1]
        l;           % pipe length (m)                                [1,1]
        nTube;       % number of tubes (-)                            [1,1]
        npass;       % number of passes (-)
        K=0;         % absolute roughness (m)                         [1,1]
        nCells;      % number of cells (-)                            [1,1]
        T1left;      % temperature left side (K)                      [1,n]
        T1right;     % temperature right side (K)                     [1,n]
        T1           % average temperature (K)                        [1,n]
        T3left;      % temperature left side (K)                      [1,n]
        T3right;     % temperature right side (K)                     [1,n]
        T3           % average temperature (K)                        [1,n]
        rho1;        % density (kg/m^3)                               [1,n]
        p1           % average pressure (Pa)                          [1,n]
        delta_p1     % pressuer loss (Pa)                             [1,n]
        p1in         % initial pressure (Pa)                          [1,1]
        h1left       % specific enthalpy left side (J/kg)             [1,n]
        h1right      % specific enthalpy right side (J/kg)            [1,n]
        h3left       % specific enthalpy left side (J/kg)             [1,n]
        h3right      % specific enthalpy right side (J/kg)            [1,n]
        Qdot1;       % heat flow (W)                                [1,n+1]
        Qdot_1       % heat flow process medium side(W)               [1,1]
        Qdot_3       % heat flow sand side(W)                         [1,1]
        QLoss_Rec    % heat flow loss (W)                             [1,n]
        alpha_i      % heat transfer coefficient inside (W/m^2*K)     [1,1]
        alpha_a = 0  % heat transfer coefficient outside (W/m^2*K)    [1,n]
        k            % overall heat transfer coefficient(W/m^2*K)     [1,n]
        F=1;         % flow characteristics crosscurrent/countercurrent [1,1]                                                                                                          
        lambda_st=50;% heat conductivity for steel tube(W/m*K)        [1,1]
        charge       % logical value (-)                              [1,n]
        QSum;
        kA;
        Atot;
        
        rec_object; % Recuperator
        comp_object; % compressor
        % sand bed
        rho_p = 2650;
        p_bed = 10^5;
        eps_por = 0.5;
        L_bed;
        w_channel;
        T_flu=0;
        A_Bett;
        h_Bett = 3;
        d_p;
        
    end
    
    methods
        function obj = hex(al_a,prop,mDot1,mDot3,di,e,l,nTube,K)
            if nargin==0
                return
            elseif nargin == 1
                obj.alpha_a = al_a;
            else
                obj.prop = prop;
                obj.mDot1 = mDot1;
                obj.mDot3 = mDot3;
                obj.di = di;
                obj.e = e;
                obj.l = l;
                obj.nTube = nTube;
                obj.K = K;
            end
        end
        
        %% TQ DIAGRAM
        % plot TQ-diagram

        function TQ_diagram(obj)
            if obj(1).process_medium == "water"
                obj(1).prop = IF97();
            elseif obj(1).process_medium == "air"
                obj(1).prop = DryAir();
            elseif obj(1).process_medium == "sCo2"
                obj(1).prop = CO2();
            end

            spec = obj(1).prop;
            len = sum([obj.nCells]);
            Qdot = [obj.Qdot_start];
            T1in = [obj.T1_start];
            a = T1in ~= 0;
            T1in = T1in(a);
            p = [obj.p_start];
            p = p(a);
            rho1in = 1./spec.v(p, T1in);
            T3in = [obj.T3_start];
            b = T3in ~= 0;
            T3in = T3in(b);

            if T1in > T3in  %charge 
                [obj.charge]=deal(true); 
                Qdot = [0, Qdot];
            else           %discharge
                [obj.charge]=deal(false);
                Qdot = [0, -Qdot];
   
            end
            step = 0;
            for n = 1:size(obj,2)
                QDot3(step+1:step+obj(n).nCells+1)=linspace((Qdot(n)),(Qdot(n))+(Qdot(n+1)),obj(n).nCells+1);
                mdot3(step+1:step+obj(n).nCells) = obj(n).mDot3;
                step = step + obj(n).nCells;
            end

            H3In=SiO2.h(T3in); % enthalpy
            if obj(1).charge
                H3r=H3In+fliplr(QDot3(1:end-1))./mdot3;
                H3l=H3In+fliplr(QDot3(2:end))./mdot3;
            else
                H3l=H3In+(QDot3(1:end-1))./mdot3;
                H3r=H3In+(QDot3(2:end))./mdot3;
            end

            T3l=SiO2.T_h(H3l); % temperature
            T3r=SiO2.T_h(H3r);
            T3m=(T3l + T3r)./2;


            P1 = NaN(1,len);
            H1l = NaN(1,len);
            H1r = NaN(1,len);
            H1In=spec.h_rhoT(rho1in,T1in); % enthalpy
            
            if obj(1).charge
                loop = 1:size(obj,2);
                step = 0;
                t=1;
            else
                loop = size(obj,2):-1:1;
                step = obj(1).nCells;
                t=-1;
            end
            for n = loop
                pos = (step+1:step+obj(n).nCells);
                QDot1m=linspace(0,(Qdot(n+1)),obj(n).nCells+1);
                if ~obj(1).charge
                    QDot1m = fliplr(QDot1m);
                end
                % water-side    
                P1in=spec.p(rho1in,T1in);
                P1(pos)=repmat(P1in,1,obj(n).nCells); % pressure
                H1l(pos)=H1In-QDot1m(1:end-1)./obj(n).mDot1;
                H1r(pos)=H1In-QDot1m(2:end)./obj(n).mDot1;
                T1l=spec.T_ph(P1,H1l); % temperature
                T1r=spec.T_ph(P1,H1r);
                if obj(1).charge
                    H1In = H1r(step+obj(n).nCells);
                else
                    H1In = H1l(step+1);
                end
                step = step + obj(n).nCells*t;
            end
            P1 = P1(~isnan(P1));
            H1l = H1l(~isnan(H1l));
            H1r = H1r(~isnan(H1r));
            T1l = T1l(~isnan(T1l));
            T1r = T1r(~isnan(T1r));

            T1m=(T1l+T1r)./2;
            H1m=(H1l+H1r)./2;

            if obj(1).process_medium == "water"
                Rho1 = hex.rho(P1,H1m); %für H2O
            else
                Rho1 = 1./spec.v(P1,T1m); % für CO2
            end

           % plot 
            if T1in > T3in % charge
                p_1 = plot(abs(QDot3),[[T1l-273.15,T1r(end)-273.15];[T3l(1)-273.15,T3r-273.15]]);
            else            % discharge
                p_1 = plot(abs(QDot3),[([T1l(1)-273.15,T1r-273.15]);([T3l-273.15,T3r(end)-273.15])]);
            end

            if obj(1).process_medium == "water"
                legend(p_1 ,{"water", 'sand'})
            elseif obj(1).process_medium == "air"
                legend(p_1,{"water", 'air'})
            elseif obj(1).process_medium == "CO2"
                legend(p_1,{"water", 'CO2'})
            end
            xlabel("QDot in W")
            ylabel("T in °C")

            step = 0;
            for n = 1:size(obj,2)
                obj(n).h1left = H1l(step+1:step+obj(n).nCells);
                obj(n).h1right = H1r(step+1:step+obj(n).nCells);
                obj(n).T1left = T1l(step+1:step+obj(n).nCells);
                obj(n).T1right = T1r(step+1:step+obj(n).nCells);
                obj(n).T1 = T1m(step+1:step+obj(n).nCells);
                obj(n).h3left = H3l(step+1:step+obj(n).nCells);
                obj(n).h3right = H3r(step+1:step+obj(n).nCells);
                obj(n).T3left = T3l(step+1:step+obj(n).nCells);
                obj(n).T3right = T3r(step+1:step+obj(n).nCells);
                obj(n).T3 = T3m(step+1:step+obj(n).nCells);
                obj(n).rho1 = Rho1(step+1:step+obj(n).nCells);
                obj(n).Qdot1 = QDot3(step+1:step+obj(n).nCells+1);
                obj(n).p1 = P1(step+1:step+obj(n).nCells);
                obj(n).p1in = P1in;
                step = step + obj(n).nCells;
            end
        end

        %% HEAT EXCHANGER CALCULATION
        % DESCRIPTIVE TEXT 
        function hex_calculation(obj,nRec,idx)
            spec = obj(1).prop;
            % get values
            len = sum([obj.nCells]);
            Rho1 = [obj.rho1];
            P1 = [obj.p1];
            T1m = [obj.T1];
            T3m = [obj.T3];
            if obj(1).charge
                T1in=[obj.T1left];
                T3in=[obj.T3right];
                h1in=[obj.h1right];
                h3in=[obj.h3left];
            else
                T1in=[obj.T1right];
                T3in=[obj.T3left];
                h1in=[obj.h1right];
                h3in=[obj.h3left];
            end

            if obj(1).flow_type == "crosscurrent"
                obj(1).F = 0.869;
            else
                obj(1).F = 1;
            end
           
            step = 0;
            for i = 1:size(obj,2)
                L(step+1:step+obj(i).nCells) = obj(i).l;
                D(step+1:step+obj(i).nCells) = obj(i).di;
                Lamb(step+1:step+obj(i).nCells) = obj(i).lambda_st;
                mdot1(step+1:step+obj(i).nCells) = obj(i).mDot1; 
                mdot3(step+1:step+obj(i).nCells) = obj(i).mDot3;
                N(step+1:step+obj(i).nCells) = obj(i).nTube;
                K_pipe(step+1:step+obj(i).nCells) = obj(i).K;  
                x(step+1:step+obj(i).nCells) = linspace(obj(i).l./obj(i).nCells,obj(i).l,obj(i).nCells); 
                QDot1(step+1:step+1+obj(i).nCells) = obj(i).Qdot1;
                Da(step+1:step+obj(i).nCells) = obj(i).di + 2*obj(i).e;
                P1In(step+1:step+obj(i).nCells) = obj(i).p_start;
                lcell(step+1:step+obj(i).nCells)=obj(i).l/obj(i).nCells;
                Ncells(step+1:step+obj(i).nCells) = obj(i).nCells;
                F1(step+1:step+obj(i).nCells) = obj(i).F;
                step = step + obj(i).nCells;
            end
            x(x<D) = D(x<D); 

            A1=Da.*pi.*L./Ncells.*N;  % cross-sectional area over all pipes
  
            if ~obj(1).charge %discharge
                T1in=fliplr(T1in);
                T3in=fliplr(T3in);
                h1in=fliplr(h1in);
                h3in=fliplr(h3in);      
                T1m=fliplr(T1m);
                T3m=fliplr(T3m);
                P1=fliplr(P1);
                Rho1=fliplr(Rho1);
                mdot1=fliplr(mdot1);
            end

            % h_bed = obj(1).h_Bett; %m
            % dp_bed = FluBed.deltaP(h_bed,obj(1).eps_por,obj(1).rho_p);
            % % p_bed = 10^5+dp_bed*0.7;
            % 
            % 
            % dp_G = dp_bed*1.7;
            % p1_G = 10^5;
            % p2_G = 10^5+ dp_G;
            % T1_G = 25 + 273.15;
            % eta_G = 0.7;
            % Tc_in = T1_G+eta_G.^-1.*T1_G.*((p2_G/p1_G).^(0.4/1.4)-1);

            
            %compressor
            if nRec ~= 1000
                comp = compressor();
                comp.h_bed = obj(1).h_Bett;
                ABett = sum([obj.A_Bett])/nRec;
                comp.A_bed = ABett;
                comp.dp = obj(1).d_p;
                comp.rho_p = obj(1).rho_p;
                comp.calculate(1.7) %facor: dp = 1.7*dp_bed
                obj(1).comp_object = comp;
            end
           
            err=repmat(10,1,5);
            QDot=0;
            count=0;
            isConv=false;
            disp([obj.name]);
           
            while err(end)>1 || ~isConv

                %heat transfer coefficient
                al_i = alpha_in(); 
                if mean(al_i) > 50000
                    al_i = repmat(50000,1,len);
                end
                if obj(1).alpha_a == 0
                    al_a = alpha_out();
                else
                    al_a = obj(1).alpha_a;
                end

                % al_a = 800;
                % al_a = 600;
                k1 = ((1./al_a) + (Da./(2.*Lamb)).*log(Da./D)+Da./(al_i.*D)).^-1;  
                dT = (T1m - T3m);

                %Wärmestrom pro Zelle
                QDoti=F1.*k1.*A1.*dT; %Charge / Discharge 

                % NTU
                % if nRec == 100
                %     Q4 = obj.QLoss_Rec;
                %     nRec = 0;
                % end

                %Recuperator
                if nRec == 0
                    Q4 = zeros(1,len);
                elseif nRec == 1000
                    Q4 = obj.QLoss_Rec;
                else
                    step_size = fix(len/nRec);
                    rest = rem(len,nRec);
                    step = 0:step_size:len-rest;
                    step(end) = step(end) + rest;
                    if nRec == 1
                        rec1 = recuperator();
                        rec = rec1;
                    elseif nRec == 2
                        rec1 = recuperator();
                        rec2 = recuperator();
                        rec = [rec1 rec2];
                    elseif nRec == 3
                        rec1 = recuperator();
                        rec2 = recuperator();
                        rec3 = recuperator();
                        rec = [rec1 rec2 rec3];
                    elseif nRec == 4
                        rec1 = recuperator();
                        rec2 = recuperator();
                        rec3 = recuperator();
                        rec4 = recuperator();
                        rec = [rec1 rec2 rec3 rec4];
                    elseif nRec == 5
                        rec1 = recuperator();
                        rec2 = recuperator();
                        rec3 = recuperator();
                        rec4 = recuperator();
                        rec5 = recuperator();
                        rec = [rec1 rec2 rec3 rec4 rec5];
                    end

                    rec_d = recuperator();
                    pSand = NaN(1,len);
                    pSand(:) = comp.p_bed;
                    for r=1:nRec
                        Th_in=mean(T3m(step(r)+1:step(r+1)));
                        mDot4_Rec = hex.mflu(comp.dp,comp.rho_p,comp.p_bed,Th_in,ABett,4);
                        rec_d.Th_in = Th_in;
                        rec_d.Tc_in = comp.Tout;
                        rec_d.mDot4 = mDot4_Rec;
                        rec_d.desing(comp.pout);
                        rec(r).recalculate(rec_d,mDot4_Rec,Th_in,comp.Tout,comp.pout);
                        ha1 = DryAir.h(comp.pout,rec(r).Tc_out);
                        ha2 = DryAir.h(pSand(step(r)+1:step(r+1)),T3m(step(r)+1:step(r+1)));
                        Q4(step(r)+1:step(r+1)) = rec(r).mDot4*(ha2-ha1)./(step(r+1)-step(r));
                    end 
                    obj(1).rec_object = rec;
                end
                QDot3 = QDoti - Q4;

                %energy balance 
                h1out=h1in-QDoti./mdot1;
                h1in(2:end)=h1out(1:end-1);
                
                T1out=spec.T_ph(P1,h1out); 
                
                %Consistency check
                if obj(1).charge
                    check=T1out<T3in;
                else
                    check=T1out>T3in;    
                end
                T1out(check)=T3in(check);
                
                T1in(2:end)=T1out(1:end-1);
                T1m=mean([T1in;T1out],1);
               
                h3out=h3in+QDot3./mdot3;
                h3in(1:end-1)=h3out(2:end);
                
                T3out=SiO2.T_h(h3out);
                
                %Consistency check
                if obj(1).charge
                    check=T3out>T1in;
                else
                    check=T3out<T1in;    
                end
                T3out(check)=T3out(check);
                
                T3in(1:end-1)=T3out(2:end);
                T3m=mean([T3in;T3out],1);

                % pressure loss
                delta = deltap_i();
                P1 = P1In -  delta;

                if obj(1).process_medium == "water"
                    Rho1 = hex.rho(P1,mean([h1in;h1out])); %für H2O
                else
                    Rho1 = 1./spec.v(P1,T1m); % für CO2
                end
              
                % termination condition
                QDotsum=sum(QDot3);
              
                err=circshift(err,-1);
                err(end)=abs(QDotsum-QDot);
                QDot=QDotsum;

           
                % output information
                if mod(count,10)==0


                    switch idx
                        case 1
                            disp(['count=',num2str(count),', err=',num2str(err(end)),' QDot=',num2str(QDot)]);
                            
                        case 2
                            disp(['count=',num2str(count),', err=',num2str(err(end)),' TH2Oout=',num2str(obj.T1(end))]);    
                    end
                end
                
                % termination if 10000 iteration
                count=count+1;
                if count>10000
                    break
                end

                %convergence condition
                derr=diff(err);
                isConv=(all(derr<0) && err(end)/err(1)>0.7) || max(abs(derr))<0.2;
            end 
            % end of iteration
            %% save values
            
            if ~obj(1).charge
                T1m = fliplr(T1m);
                T3m = fliplr(T3m);
                T1in = fliplr(T1in);
                T1out = fliplr(T1out);
                T3out = fliplr(T3out);
                T3in = fliplr(T3in);
                h1in = fliplr(h1in);
                h1out = fliplr(h1out);
                h3out = fliplr(h3out);
                h3in = fliplr(h3in);
                P1 = fliplr(P1);
                delta = fliplr(delta);
                QDoti = fliplr(QDoti);
                QDot3 = fliplr(QDot3);
                % Q4 = fliplr(Q4);
                % obj.deltapi_ausgabe = fliplr(obj.deltapi_ausgabe);
                % obj.QLoss_Rec = Q4;


            end
            step = 0;
            
            % obj.k = k1;
            
            for t = 1:size(obj,2)
                if obj(1).charge
                    obj(t).T1 = T1m(step+1:step+obj(t).nCells);
                    obj(t).T3 = T3m(step+1:step+obj(t).nCells);

                    obj(t).T1left = T1in(step+1:step+obj(t).nCells);
                    obj(t).T1right = T1out(step+1:step+obj(t).nCells);
                    obj(t).T3left = T3out(step+1:step+obj(t).nCells);
                    obj(t).T3right = T3in(step+1:step+obj(t).nCells);

                    obj(t).h1left = h1in(step+1:step+obj(t).nCells);
                    obj(t).h1right = h1out(step+1:step+obj(t).nCells);
                    obj(t).h3left = h3out(step+1:step+obj(t).nCells);
                    obj(t).h3right = h3in(step+1:step+obj(t).nCells);
                else
                    obj(t).T1 = T1m(step+1:step+obj(t).nCells);
                    obj(t).T3 = T3m(step+1:step+obj(t).nCells);

                    obj(t).T1left = T1out(step+1:step+obj(t).nCells);
                    obj(t).T1right = T1in(step+1:step+obj(t).nCells);
                    obj(t).T3left = T3in(step+1:step+obj(t).nCells);
                    obj(t).T3right = T3out(step+1:step+obj(t).nCells);

                    obj(t).h1left = h1out(step+1:step+obj(t).nCells);
                    obj(t).h1right = h1in(step+1:step+obj(t).nCells);
                    obj(t).h3left = h3in(step+1:step+obj(t).nCells);
                    obj(t).h3right = h3out(step+1:step+obj(t).nCells);
                end

                obj(t).p1 = P1(step+1:step+obj(t).nCells);
                obj(t).delta_p1 = delta(step+1:step+obj(t).nCells);
                obj(t).delta_p1 = delta(step+1:step+obj(t).nCells);
                obj(t).Qdot_1 = QDoti(step+1:step+obj(t).nCells);
                obj(t).Qdot_3 = QDot3(step+1:step+obj(t).nCells);
                obj(t).alpha_i = mean(al_i(step+1:step+obj(t).nCells));
                obj(t).alpha_a = al_a;
                obj(t).k = mean(k1(step+1:step+obj(t).nCells));
                obj(t).QLoss_Rec = Q4(step+1:step+obj(t).nCells);
                obj(t).QSum = sum(QDoti(step+1:step+obj(t).nCells));

                step = step + obj(t).nCells;
            end
                
            %% HEAT TRANSFER COEFFICIENT INSIDE
            % assumption: 
            % single phase: constant heat flow density
            % two phase: vertical pipes
        
            function alpha_1=alpha_in()
                persistent g 
                if isempty(g)
                    g = 9.80665; %acceleration of gravity
                end
                is2phase = spec.is2phase(Rho1,T1m);
                % value = zeros(1,len);
                % is2phase = value == 1;
                alphainside = NaN(1,len); 

                a = find(is2phase==0);

                A = D.^2.*pi/4;
                mDot = mdot1./(A.*N);

                if a > 0
                    eta = spec.my(Rho1,T1m); % viscosity
                    lambda1 = spec.lambda(Rho1,T1m); %heat conductivity
                    u = mDot./Rho1; %velocity
                    Pr = spec.Pr(Rho1,T1m); % Pr-number
                    Re = u.*D.*Rho1./eta; % Re-number
                
                    %%% SINGL PHASE
                    % G1-3 Wärmeübertragung bei laminarer Strömung durch Rohre (S.787)
                    % G1-4 Wärmeübertragung bei turbulenter Strömung durch Rohre (S.788)

                    Nu_x = NaN(1, len); %local Nu-number
                    lam = Re < 2300; % laminar flow
                    Nu_x2 = 1.302 .* (Re.*Pr.*D./x).^(1/3);
                    Nu_x(lam) = (4.364.^3 + 1 + (Nu_x2(lam) - 1).^3).^(1./3);
                    tur = Re > 10^4; % turulent flow
                    Xi = (1.8 .* log10(Re)-1.5).^(-2);
                    Nu_x(tur) = (Xi(tur)./8).*(Re(tur)-10^3).*Pr(tur)./(1+12.7.*(Xi(tur)./8).^0.5.*(Pr(tur)-1)).*(1+1/3.*(D(tur)./x(tur)).^(2/3)); %neue Formel

                    tra = 2300 < Re & Re < 10^4; %transition region
                    Nu_x2L = 1.302 .* (2300.*Pr(tra).*D(tra)./x(tra)).^(1./3);
                    Nu_x3L = 0.462*Pr(tra).^(1/3).*(2300.*D(tra)./x(tra)).^(1/2);
                    Nu_xL = (4.364^3+(Nu_x2L-0.6).^3+(Nu_x3L).^3).^(1/3);
                    Nu_xT = (0.0308./8).*10^4.*Pr(tra)./(1+12.7.*(0.0308./8).^(1/2).*(Pr(tra).^(2./3)-1)).*(1+(1/3).*(D(tra)./x(tra)).^(2/3));
                    Gamma = (Re(tra) - 2300)./(10.^4 -2300);
                    Nu_x(tra) = (1-Gamma) .* Nu_xL + Gamma .* Nu_xT;

                    alphainside(~is2phase) = Nu_x(~is2phase) .* lambda1(~is2phase) ./ D(~is2phase); % heat transfer coefficient
                end

                % H3.4.2 Blasensieden reiner Stoffe in durchströmten Rohren (S.919)
                persistent cf q0 d0 Ra0 pc alpha_0                    
                if isempty(cf)
                    cf = spec.cf;    %properties of the fluid                                                    
                    q0 = spec.q0;  %substance-specific values
                    d0 = spec.d0;
                    Ra0 = spec.Ra0;
                    pc = spec.p_c;
                    alpha_0 = spec.alpha_0;
                end
                Ra_in = 0.5 * 10^-6; %Angenommen
                p_st = P1./pc;
                n = 0.8 - 0.1.*10.^(0.76 .* p_st);

                q = abs(QDot1(2:end)./(D.^2*pi./4));

                A = NaN(1,len);
                B = NaN(1,len);
                C = NaN(1,len);

                A(is2phase) = cf.*(q(is2phase)./q0).^n(is2phase);
                B(is2phase) = 2.816 .* p_st(is2phase) .^0.45 + (3.4 + 1.7./(1-p_st(is2phase).^7)).*p_st(is2phase).^3.7;
                C(is2phase) = (d0./D(is2phase)).^0.4 .* (Ra_in./Ra0).^0.133;

                alphainside(is2phase) = A(is2phase).*B(is2phase).*C(is2phase).*alpha_0; % heat transfer coefficient 
                alpha_1 = alphainside;
            end

            %% HEAT TRANSFER COEFFICIENT OUTSIDE
            function alpha_3 = alpha_out()
                persistent lambda40 baseAlpha
                if isempty(lambda40)
                    lambda40=tppAir('lambda',40+273.15)^0.6;
    
                    d_sv=120e-6;
                    baseAlpha=1200.*d_sv.^-0.36./150e-6.^-0.36;
                end
                alpha=baseAlpha.*tppAir('lambda',T3m).^0.6./lambda40;
                alpha_3 = alpha;
            end   
            
            %% PRESSURE LOSS
            function deltap = deltap_i()
                persistent g 
                if isempty(g)
                    g = 9.80665; %acceleration of gravity
                end

                is2phase = spec.is2phase(Rho1,T1m);
                % value = zeros(1,len);
                % is2phase = value == 1;

                rho_l = NaN(1,len);
                rho_g = NaN(1,len);
                xDot = NaN(1,len);
                eta_l = NaN(1,len);
                eta_g = NaN(1,len);
                eta = NaN(1,len);
    
                rho_l(is2phase) = spec.rho_satL(T1m(is2phase));
                rho_g(is2phase) = spec.rho_satV(T1m(is2phase));
                xDot(is2phase) = spec.x(Rho1(is2phase),T1m(is2phase));
                eta_l(is2phase) = spec.my(rho_l(is2phase),T1m(is2phase));
                eta_g(is2phase) = spec.my(rho_g(is2phase),T1m(is2phase));
                eta(~is2phase) = spec.my(Rho1(~is2phase),T1m(~is2phase)); 
              
                A = D.^2.*pi/4; 
                mDot = mdot1./(N.*A); %mass flow density
                u = mDot./Rho1; %velocity
                Re = u(~is2phase).*Rho1(~is2phase).*D(~is2phase)./eta(~is2phase);         
                deltap_i = NaN(1,len);
        
                %%% SINGLE PHASE
                % Druckverlust in durchströmten Rohren 
                % 12.Auflage
                f = NaN(1, size(Re,2)); %empty friction factor  array
                lam=Re<3000; % laminar flow
                f(lam)=64./Re(lam); %eq (4)
                tur = Re>=3000; %turbulent flow              
                fP = 64 ./ Re(tur);% start value
                d = 1;
                K_sp = K_pipe(~is2phase);
                D_sp = D(~is2phase);
                while d > 0.000001 % eq (10)
                    z = K_sp(tur)./D_sp(tur)./3.71+ 2.51./(Re(tur).*fP.^0.5);
                    n = -2 .* log10(z);
                    f_new = (1./n).^2;
                    d = sum(abs(fP - f_new))/size(Re,2);
                    fP = f_new;
                end
                f(tur) = fP;
                deltap_i(~is2phase) = f .* (lcell(~is2phase)./D(~is2phase))...
                    .* (Rho1(~is2phase) .* u(~is2phase).^2 ./ 2); % eq (1)

                %%% TWO PHASE VAPORIZATION
                % Strömungssieden – Druckverlust in durchströmten Verdampferrohren
                % VDI 12.Auflage 
                % This chapter has been revised in the new VDI Heat Atlas
                % edition. Therefore new equations are used.
                
                % friction pressure loss
                % heterogenes Modell nach Friedel
                dh = D; %hydraulic diameter               
                sigma = spec.sigma(T1m); %surface tension
                Re_LO = mDot.*dh./(eta_l); 
                Re_GO = mDot.*dh./(eta_g);
                ReLO_s = Re_LO <= 1055;
                ReLO_b = Re_LO > 1055;
                ReGO_s = Re_GO <= 1055;
                ReGO_b = Re_GO > 1055;

                Zeta_LO = NaN(1,len);
                Zeta_GO = NaN(1,len);

                Zeta_LO(ReLO_s) = 64 .* Re_LO(ReLO_s) .^-1;
                Zeta_GO(ReGO_s) = 64 .* Re_GO(ReGO_s) .^-1;
                Zeta_LO(ReLO_b) = 0.86859.*log((Re_LO(ReLO_b)./(1.964.*log(Re_LO(ReLO_b))-3.8215))).^-2;
                Zeta_GO(ReGO_b) = 0.86859.*log((Re_GO(ReGO_b)./(1.964.*log(Re_GO(ReGO_b))-3.8215))).^-2;

                Fr_LO = mDot.^2.*(g.*dh.*rho_l.^2).^-1;
                We_LO = mDot.^2.*dh.*(rho_l.*sigma).^-1;
                
                % two-phase multiplier 
                % eq (10) -(11)
                Ae = (1-xDot).^2 + xDot.^2.*(rho_l.*Zeta_GO.*(rho_g.*Zeta_LO).^-1);
                Phi_LO =  Ae + 3.43.*xDot.^0.685.*(1-xDot).^0.24.*(rho_l.*rho_g.^-1).^0.8...
                    .*(eta_g.*eta_l.^-1).^0.22.*(1-eta_g.*eta_l.^-1).^0.89...
                    .* Fr_LO.^-0.047.*We_LO.^-0.087;
                
                dp_dl = Zeta_LO.*mDot.^2.*(2.*dh.*rho_l).^-1; %gradient liquid phase
                grad_friction = Phi_LO(is2phase).*dp_dl(is2phase);
                dp_friction = grad_friction .* lcell(is2phase); %friction pressure loss

                % static pressure loss
                % Zero at horizontal pipes

                %accelaration pressure loss
                % eq (20)
                eps = xDot./rho_g .* ((1 + 0.12 .*(1-xDot)) .* (xDot./rho_g + (1-xDot)./rho_l)...
                    + 1.18.*(1-xDot).*(g.*sigma.*(rho_l-rho_g)).^0.25./(mDot.*rho_l.^0.5)).^-1; %gas volume fraction
                int = xDot.^2./(eps.*rho_g) + (1-xDot).^2./((1-eps).*rho_l);
                dif = [(int(2:end) - int(1:end-1)) 0];
                dp_ac = mDot.^2.*lcell.*dif;
                dp_ac(isnan(dp_ac))=0;

                %%% ADDITIONAL PRESSURE LOSS
                %FDBR-Handbuch März 2017 
                % Kapitel 9

                % Cross section narrowing
                Zeta_s = 0.05; %(9-32)

                % Cross section expansion
                D2 = 40 * 10^-3;
                Zeta_SD = (1-D./D2).^2;

                
                %Pipe bends
                del = 180; %180° bend
                rKR = D;
                xKR = rKR./D;
                c1 = 9.3.*exp(-0.06.*xKR);
                c2 = 0.0788.*tanh(0.8.*xKR) + 0.00124.*xKR;
                c3 = 15000.*exp(-2.7.*xKR) + 1780 .* exp(-0.0234.*(xKR-8).^2);
                c4 = 4.95 + 4.042.*exp(-0.01.*abs(xKR-7.5).^3);
                fKR = c1 .* tanh(c2./c2.*(del+10.05)) + c3.*(del./100).^4.*exp(-c4.*del./100);
                if sum(is2phase) ~= 0
                    f = 0.025;
                end
                Zeta_KR = f .* fKR;

                % dp_add = zeros(1,Ncells(1)+Ncells(end));
                dp_add = zeros(1,len);
                dp_add(1) = Zeta_s(1).*Rho1(1)./2.*u(1).^2;
                dp_add(end) = Zeta_SD(end) .* Rho1(end)./2.*u(end).^2;
       
                %180° Pipebends
                bend = round(len./sum([obj.npass]))-1; %ACHTUNG!!
                for n=1:sum([obj.npass])-1
                    dp_add(bend*n) = Zeta_KR(bend*n) .* Rho1(bend*n)./2.*u(bend*n).^2;
                end

                % whole loss of pressure
                % eq (1)
                deltap_i(is2phase) = dp_friction + dp_ac(is2phase);
                deltap_i = deltap_i + dp_add;
                deltap = cumsum(deltap_i);
                
            end 
        end

        function plot_Tx(obj,type)
            if nargin < 2
                type = 1;
            end
            T1l = [obj.T1left];
            T1r = [obj.T1right];
            T3l = [obj.T3left];
            T3r = [obj.T3right];
            step = 0;
            step2 = 0;
            XI = NaN(1,sum([obj.nCells]));
            for n = 1:size(obj,2)
                XI(step+1:step+obj(n).nCells) = step2 + linspace(obj(n).L_bed./obj(n).nCells,obj(n).L_bed,obj(n).nCells);
                if type == 2
                    p_1 = plot([step2 XI(step+1:step+obj(n).nCells)], [T1l(step+1) T1r(step+1:step+obj(n).nCells)]-273.15,'Color',[0 0.4470 0.7410]);
                    hold on
                end
                step = step + obj(n).nCells;
                step2 = XI(step);
            end
            XII = [0 XI];
            if type == 1
                p_1 = plot(XII,[T1l(1) T1r]-273.15,'Color',[0 0.4470 0.7410]);
                hold on
            end
            p_3 = plot(XII,[T3l(1) T3r]-273.15,'Color',[0.9290 0.6940 0.1250]);
            hold off
            xlabel("x in m")
            ylabel("T in °C")
            if obj(1).process_medium == "water"
                legend([p_1 p_3],{"water", 'sand'})
            elseif obj(1).process_medium == "air"
                legend([p_1 p_3],{"water", 'air'})
            elseif obj(1).process_medium == "CO2"
                legend([p_1 p_3],{"water", 'CO2'})
            end
        end

        function plot_TQ(obj, type)
            if nargin < 1
                type = 1;
            end
            Q = cumsum(abs([obj.Qdot_1]*10^-6));
            T1l = [obj.T1left];
            T1r = [obj.T1right];
            T3l = [obj.T3left];
            T3r = [obj.T3right];   
            QII = [0 Q];

            step = 0;
            step2 = 0;
            % XII = NaN(1,sum([Gaston.nCells]));
            if type == 2
                for n = 1:size(obj,2)
                    QI = Q(step+1:step+obj(n).nCells);
                    p_1 = plot([step2 QI],[T1l(step+1) T1r(step+1:step+obj(n).nCells)]-273.15,'Color',[0 0.4470 0.7410]);
                    hold on
                    step = step + obj(n).nCells;
                    step2 = QI(end);
                end
            elseif type == 1
                p_1 = plot(QII,[T1l(1) T1r]-273.15,'Color',[0 0.4470 0.7410]);
                hold on
            end

            p_3 = plot(QII,[T3l(1) T3r]-273.15,'Color',[0.9290 0.6940 0.1250]);
            hold off
            xlabel("Q in MW")
            ylabel("T in °C")
            if obj(1).process_medium == "water"
                legend([p_1 p_3],{"water", 'sand'})
            elseif obj(1).process_medium == "air"
                legend([p_1 p_3],{"water", 'air'})
            elseif obj(1).process_medium == "CO2"
                legend([p_1 p_3],{"water", 'CO2'})
            end
            
        end
    end
    
    % density characteristic
    methods (Static)
        function rho=rho(p,h)
            spec = IF97;
            hsatl=spec.h(p,NaN,0);
            hsatv=spec.h(p,NaN,1);
            is2phase=hsatl<=h & h<=hsatv;
            
            x=NaN(size(h));
            x(is2phase)=(h(is2phase)-hsatl(is2phase))./(hsatv(is2phase)-hsatl(is2phase));
            
            T=NaN(size(h));
            T(~is2phase)=spec.T_ph(p(~is2phase),h(~is2phase));
            
            rho=1./spec.v(p,T,x);

        end

        function mf = mflu(dp,rhop,p,T,A,flu)
            uf = FluBed.wmf(dp,rhop,p,T);
            rho_air = DryAir.rho(p,T);
            mf = rho_air*uf*A*flu;
        end

    end
 end




