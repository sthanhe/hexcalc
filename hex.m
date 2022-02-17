classdef hex < handle
    %RECALCULATION HEAT EXCAHNGER
    %references: VDI Wärmeatlas, 11.Auflage, SpringerVieweg
                  
     % charge/discharge
     %1: H2O-Side
     %3: Sand-Side
     
    properties
        mDot1;       % mass flow (kg/s)                               [1,1]
        mDot3;       % mass flow (kg/s)                               [1,1]
        di;          % inner diameter (m)                             [1,1]
        e;           % wall thickness (m)                             [1,1]
        l;           % pipe length (m)                                [1,1]
        N;           % number of pipes (-)                            [1,1]
        Ra;          % mean roughness index (m)                       [1,1] 
        K=0;         % absolute roughness (m)                         [1,1]
        nCells=100;  % number of cells (-)                            [1,1]
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
        h1In         % initial specific enthalpy (J/kg)               [1,1]
        h3In         % initial specific enthalpy (J/kg)               [1,1]
        h1left       % specific enthalpy left side (J/kg)             [1,n]
        h1right      % specific enthalpy right side (J/kg)            [1,n]
        h3left       % specific enthalpy left side (J/kg)             [1,n]
        h3right      % specific enthalpy right side (J/kg)            [1,n]
        Qdot1;       % heat flow (W)                                [1,n+1]
        Qdot         % heat flow (W)                                  [1,1]
        alpha_i      % heat transfer coefficient inside (W/m^2*K)     [1,n]
        alpha_a      % heat transfer coefficient outside (W/m^2*K)    [1,n]
        k            % overall heat transfer coefficient(W/m^2*K)     [1,n]
        F=0.869;     % flow characteristics (-)                       [1,1]                                                                                                          
        lambda_st=50;% heat conductivity for steel tube(W/m*K)        [1,1]
        charge       % logical value (-)                              [1,n]
        count_a
        err_a  
    end
    
    methods
        function obj = hex(mDot1,mDot3,di,e,l,N,Ra,K)
            if nargin==0
                return
            end
            obj.mDot1 = mDot1;
            obj.mDot3 = mDot3;
            obj.di = di;
            obj.e = e;
            obj.l = l;
            obj.N = N;
            obj.Ra = Ra;
            obj.K = K;
        end
        
        %% TQ DIAGRAM
        % plot TQ-diagram
        % compute and save start values
        function TQ_diagram(obj,Qdot,T1in,rho1in,T3in)
            if T1in > T3in  %charge
                obj.Qdot1=linspace(0,abs(Qdot),obj.nCells+1);
                obj.charge=true;
            else            %discharge
                obj.Qdot1=linspace(-abs(Qdot),0,obj.nCells+1);
                obj.charge=false;
            end
            
            % water-side
            obj.p1=IF97.p(rho1in,T1in); 
            obj.p1in=IF97.p(rho1in,T1in); 
            
            obj.h1In=IF97.h_rhoT(rho1in,T1in); % enthalpy
            obj.h1left=obj.h1In-obj.Qdot1(1:end-1)./obj.mDot1; 
            obj.h1right=obj.h1In-obj.Qdot1(2:end)./obj.mDot1;
            h1=(obj.h1left+obj.h1right)/2;
            
            obj.T1left=IF97.T_ph(obj.p1,obj.h1left); % temperature
            obj.T1right=IF97.T_ph(obj.p1,obj.h1right);  
            obj.T1=(obj.T1left+obj.T1right)./2;
            
            obj.p1=repmat(obj.p1,1,obj.nCells); % pressure
            
            obj.rho1=hex.rho(obj.p1,h1); % density
            
            % sand-side
            obj.h3In=SiO2.h(T3in); % enthalpy
            obj.h3left=obj.h3In+fliplr(obj.Qdot1(2:end))./obj.mDot3; 
            obj.h3right=obj.h3In+fliplr(obj.Qdot1(1:end-1))./obj.mDot3;
            
            obj.T3left=SiO2.T_h(obj.h3left); % temperature
            obj.T3right=SiO2.T_h(obj.h3right);   
            obj.T3=(obj.T3left+obj.T3right)./2;
            
            % plot
            QDot=linspace(0,abs(Qdot),obj.nCells+1);                        
            if T1in > T3in % charge
                plot(QDot,[[obj.T1left,obj.T1right(end)];[obj.T3left(1),obj.T3right]]);
            else % discharge
                plot(Qdot,[[obj.T1left(1),obj.T1right];[obj.T3left,obj.T3right(end)]]);
            end
            legend('Wasser','Sand')
            xlabel("QDot in W")
            ylabel("T in K")
        end
        
        %% HEAT TRANSFER COEFFICIENT INSIDE
        % assumption: 
        % single phase: constant heat flow density
        % two phase: vertical pipes
        
        function alpha_in(obj)
            x=linspace(obj.l./obj.nCells,obj.l,obj.nCells); 
            x(x<obj.di) = obj.di;
            eta = IF97.my(obj.rho1,obj.T1); % viscosity                    
            lambda1 = IF97.lambda(obj.rho1,obj.T1); %heat conductivity            
            u = 4.*obj.mDot1./(obj.rho1 .* obj.N .* obj.di.^2 .* pi); %velocity
            Pr = IF97.Pr(obj.rho1,obj.T1); % Pr-number
            Re = u.*obj.di.*obj.rho1./eta; % Re-number   
            is2phase = IF97.is2phase(obj.rho1,obj.T1);
            alpha = NaN(size(Re));  
            
            %%% SINGL PHASE
            % G1-3 Wärmeübertragung bei laminarer Strömung durch Rohre (S.785)
            % G1-4 Wärmeübertragung bei turbulenter Strömung durch Rohre (S.788)
            Re_si = Re(~is2phase);
            Pr_si = Pr(~is2phase);
            
            Nu_x = NaN(1, size(Re_si,2)); %local Nu-number                 
            lam = Re_si < 2300; % laminar flow
            Nu_x2 = 1.302 .* (Re_si(lam).*Pr_si(lam).*obj.di(lam)./x(lam)).^(1/3);
            Nu_x(lam) = (4.364.^3 + 1 + (Nu_x2 - 1).^3).^(1./3);
            tur = Re_si > 10^4; % turulent flow
            Xi = (1.8 .* log10(Re_si(tur))-1.5).^(-2);
            Nu_x(tur) = (Xi./8).*Re_si(tur).*Pr_si(tur)./(1+12.7.*(Xi./8).^(1/2).*(Pr_si(tur).^(2/3)-1)).*(1+(1/3).*(obj.di./x(tur)).^(2/3));
            tra = 2300 < Re_si & Re_si < 10^4; %transition region                       
            Nu_x2L = 1.302 .* (2300.*Pr_si(tra).*obj.di(tra)./x(tra)).^(1./3);
            Nu_x3L = 0.462*Pr_si(tra).^(1/3).*(2300.*obj.di(tra)./x(tra)).^(1/2);
            Nu_xL = (4.364^3+(Nu_x2L-0.6).^3+(Nu_x3L).^3).^(1/3);
            Nu_xT = (0.0308./8).*10^4.*Pr_si(tra)./(1+12.7.*(0.0308./8).^(1/2).*(Pr_si(tra).^(2./3)-1)).*(1+(1/3).*(obj.di(tra)./x(tra)).^(2/3));
            Gamma = (Re_si(tra) - 2300)./(10.^4 -2300);
            Nu_x(tra) = (1-Gamma) .* Nu_xL(tra) + Gamma .* Nu_xT(tra);
            
            alpha(~is2phase) = Nu_x .* lambda1(~is2phase) ./ obj.di; % heat transfer coefficient
            
            %%% TWO PHASE VAPORIZATION
            % H3.4.2 Blasensieden reiner Stoffe in durchströmten Rohren (S.919)
            persistent cf q0 d0 Ra0 pc alpha_0                    
            if isempty(cf)                                                 
                cf = 0.72;    %properties of the fluid                                                    
                q0 = 150000;  %substance-specific values
                d0 = 10^-2;
                Ra0 = 10^-6;
                pc = 220.64 * 10^5;
                alpha_0 = 25580;
            end
            
            p_st = obj.p1(is2phase)./pc;
            n = 0.8 - 0.1.*10.^(0.76 .* p_st);
            
            q = abs(obj.Qdot1(2:end)./(obj.di^2*pi/4));
            A = cf.*(q(is2phase)./q0).^n;
            B = 2.816 .* p_st .^0.45 + (3.4 + 1.7./(1-p_st.^7)).*p_st.^3.7;
            C = (d0./obj.di).^0.4 .* (obj.Ra/Ra0).^0.133;
            alpha(is2phase) = A.*B.*C.*alpha_0; % heat transfer coefficient 
            
            obj.alpha_i = sum(alpha)/obj.nCells; % average heat transfer coefficient
        end
            
        %% PRESSURE LOSS
        % assumption: 
        
        function deltap_i = deltap_i(obj)
            persistent g 
            if isempty(g)
                g = 9.80665; %acceleration of gravity          
            end
                   
            eta = NaN(1,size(obj.T1,2));
            is2phase = IF97.is2phase(obj.rho1,obj.T1);
            c_g = IF97.w(obj.p1(is2phase),obj.T1(is2phase),1); 
            rho_l = IF97.rho_satL(obj.T1(is2phase));
            rho_g = IF97.rho_satV(obj.T1(is2phase));
            xDot = IF97.x(obj.rho1(is2phase),obj.T1(is2phase));
            eta_l = IF97.my(rho_l,obj.T1(is2phase));
            eta_g = IF97.my(rho_g,obj.T1(is2phase));
            eta(is2phase) = xDot.*eta_g + (1-xDot).*eta_l;
            eta(~is2phase) = IF97.my(obj.rho1(~is2phase),obj.T1(~is2phase)); 
            
            A = obj.di^2*pi/4; 
            u = obj.mDot1./(obj.rho1.*A*obj.N); 
            Re = u.*obj.rho1.*obj.di./eta;         
            deltap_i = NaN(1,size(Re,2));
    
            %%% SINGL PHASE
            % L1.2 Druckverlust in durchströmten Rohren (S.1223)
            Re_si = Re(~is2phase);
            f = NaN(1, size(Re_si,2)); %empty friction factor  array
            
            lam=Re_si<3000; % laminar flow
            f(lam)=64./Re_si(lam);
            tur_B = 3*10^3 <= Re_si & Re_si < 10^4; % turbulent flow 
            f(tur_B) = 0.3164./(Re_si(tur_B).^(1/4));
            tur_K = 10^4 <= Re_si & Re_si <= 10^6;% turbulent flow 
            f(tur_K) = (1.8 .* log10(Re_si(tur_K)) - 1.5).^(-2);
            tur_P = Re_si > 10^6; % turbulent flow
            fP = 64 ./ Re_si(tur_P);% start value
            d = 1;
            while d > 0.000001 % exactness of the interation
                x = Re_si(tur_P).*(fP.^0.5);
                n = 2 .* log10(x) - 0.8;
                f_new = (1./n).^2;
                d = sum(abs(fP - f_new))/size(Re_si(tur_P),2);
                fP = f_new;
            end
            f(tur_P) = fP;
   
            x=linspace(obj.l./obj.nCells,obj.l,obj.nCells);
            x(x<obj.di) = obj.di;
            deltap_i(~is2phase) = f .* (x(1)./obj.di) .* (obj.rho1(~is2phase) .* u(~is2phase).^2 ./ 2); %prassure loss

            %%% TWO PHASE VAPORIZATION
            % H3.2 Druckverlust in durchströmten Verdampferrohren (S.903)      
            Xi = NaN(1, size(obj.T1(is2phase),2));
            mDot = obj.mDot1/obj.N;
            
            Fr = u(is2phase).^2./(g*obj.di);
            a = xDot.*rho_l./((1-xDot).*rho_g); 
            beta = 1./a;
          
            grad_friction = NaN(1, size(eta_l,2));
            eps = NaN(1, size(rho_l,2));
            
            % friction pressure loss
            dis = a <= (12.*Fr.^0.5./((1+Fr.^0.5./7))); % dispersed
            K2 = NaN(1, size(a(dis),2));%empty friction pressure loss array
            bs = beta(dis) <= 0.4;
            K2(bs) = 1 + 0.09 .* beta(bs); 
            bb = beta(dis) > 0.4;
            K2(bb) = 1./(1 - (2.97./(beta(bb).^(2./3))+1)./(6.*(1.83./(beta(bb).^(2./3))+1).*(3.43./(beta(bb).^(2./3))+1)));
            Re_ZP = mDot.*obj.di./(eta(dis).*(1-xDot(dis).*(1-eta_g(dis)./eta_l(dis))));
            
            %friction factor
            xid = 0.1; %Startwert
            d = 1;
            while d > 0.000001
                n = -2.*log10(obj.K./(3.7.*obj.di))+2.51./(Re_ZP.*xid.^0.5);
                xi_new = (1./n).^2;
                d = sum(abs(xid - xi_new))/size(xi_new,2);
                xid = xi_new;
            end
            Xi(dis) = xid;
            
            % friction pressure loss
            A = 1 + xDot(dis).*(rho_l(dis)./rho_g(dis)-1);
            B = 1 - xDot(dis).*(rho_l(dis)./rho_g(dis) -1).*(K2 - 1);
            grad_friction(dis) = Xi(dis).*mDot.^2./(2.*rho_g(dis).*obj.di).*A.*B;
            
            %continous gas phase 
            con = a > (12.*Fr.^0.5./((1+Fr.^0.5./7)));
            
            Re_l = mDot.*obj.di.*(1-xDot(con))./(eta_l(con));
            Re_g = mDot.*obj.di.*xDot(con)./(eta_g(con));
            Fr_l = mDot.^2 .* (1-xDot(con)).^2 ./ (rho_l(con).*g.*obj.di);
            psi = (1-xDot(con))./xDot(con) .* (Re_l.*Fr_l).^(-1./6).*(rho_l(con)./rho_g(con)).^(-0.9) .* (eta_l(con)./eta_g(con)).^(-0.5);
            
            
            if obj.K/obj.di < 5*10^(-4) %smooth pipes
            eps_1 = 1.71.*psi.^0.2 .* ((1-xDot(con))./xDot(con)).^0.15 .* (rho_g(con)./rho_l(con)).^0.1;
            else %not smooth pipes
            eps_1 = 1.71.*psi.^0.2 .* ((1-xDot(con))./xDot(con)).^0.5 .* (rho_g(con)./rho_l(con)).^0.1 .* (5.*10.^4./(obj.K./obj.di)).^0.13;
            end
            eps_2 = 9.1 .* psi;
            eps_f = (eps_1.^-3 + eps_2.^-3).^(-1/3);
             
            %Barnett number
            gamma_E = (1 + 6.67./(((1-xDot(con))./xDot(con)).^0.45 .* (1+3.*xDot(con).^4).*(eta_l(con)./eta_g(con)-1).^0.25)).^(-1);
            gamma_F = 1 - (1 + (1-xDot(con)).*rho_g(con)./(xDot(con).*eps_f.*rho_l(con))).^(-1.19);
            E = 1.857 + 0.814 .* log10((mDot.*xDot(con)./(rho_g(con).*c_g(con))).^2 .* (1 + 4575 .* rho_g(con).^2./rho_l(con).^2));
            
            phi = (1./(1-(1-E).*gamma_F-E.*gamma_E)).^2;
            
            %frictionnumber
            xic = 0.1; % start value
            d = 0.1;
            while d > 0.000001
                n = 2*log10(Re_g.*xic.^0.5) - 0.8;
                Xi_new = (1./n).^2;
                d = sum(abs(xic - Xi_new))/size(Xi_new,2);
                xic = Xi_new;
            end
            Xi(con) = xic;
            grad_friction(con) = Xi(con) .* mDot.^2./(2.*rho_g(con).*obj.di).*phi;
            dp_friction = grad_friction .* x(1);
            
            %STATIC PRESSURE LOSS
            %dispersed
            K_dis = 1 - (30.4./beta(dis).^2.3 + 11)./(60.*(1.6./beta(dis).^2.3 +1).*(3.2./beta(dis).^2.3+1));
            eps(dis) = K_dis./(1+beta(dis));
            
            %continous gas phase
            H = NaN(1, size(a(con),2));
            eps_con = eps(con);
            beta_con = beta(con);
            xDot_con = xDot(con);
            eta_l_con = eta_l(con);
            eta_g_con = eta_g(con);
            rho_l_con = rho_l(con);
            rho_g_con = rho_g(con);
            
            aa = a(con) > 10^4;
            H(aa) = beta_con(aa)./(1+beta_con(aa));
            eps_con(aa) = 1-H(aa);
            ab = 10^4 > a(con) & a(con) > 500;
            H(ab) = 464 .* beta_con(ab).^(5/3);
            eps_con(ab) = 1-H(ab);
            ac = a(con) < 500 & a(con) > 12.*Fr(con).^0.5./(1+Fr(con).^0.5./7); 
            X_u = ((1-xDot_con(ac))./xDot_con(ac)).^(7./8) .* (eta_l_con(ac)./eta_g_con(ac)).^(1/8) .* (rho_g_con(ac)./rho_l_con(ac)).^(1/2);
            H_2 = X_u./(1+X_u);
            eps_con(ac) = 0.1; %start value
            while d < 0.000001
                H_1 = exp(2-0.1335.*log(eta_l_con(ac)/eta_g_con(ac)) + (1.1 - 0.08534*log(eta_l_con(ac)/eta_g_con(ac))))*log(eps_con(ac));
                H(ac) = (H_1.^-3 + H_2.^-3).^(-1/3);
                eps_new = 1 - H(ac);
                d = sum(abs(eps_con(ac) - eps_new))/size(eps_new,2);
                eps_con(ac) = eps_new;
            end 
            eps(con) = eps_con;
            grad_static = (rho_l.*(1-eps) + rho_g.*eps)*g;
            dp_static = grad_static .* x(1);
            
            %ACCELERATION PRESSURE LOSS
            dp_acceleration = NaN(1, size(dp_static,2));
            if size(eps,2) > 0
                for pos = 2:(size(dp_static,2))
                    part_1 = xDot(pos-1).^2./(eps(pos-1).*rho_g(pos-1)) + (1-xDot(pos-1)).^2./((1-eps(pos-1)).*rho_l(pos-1));
                    part_2 = xDot(pos).^2./(eps(pos).*rho_g(pos)) + (1-xDot(pos)).^2./((1-eps(pos)).*rho_l(pos));
                    dp_acceleration(pos) = mDot.^2.*(part_2 - part_1);
                end 
                dp_acceleration(1) = xDot(1).^2./(eps(1).*rho_g(1)) + (1-xDot(1)).^2./((1-eps(1)).*rho_l(1));
            else
                dp_acceleration = NaN(1,size(eps,2)); 
            end
            
            %WHOLE PRESSURE LOSS
            deltap_i(is2phase) = (dp_friction + dp_static + dp_acceleration);    
            
            for pos = 1:size(deltap_i,2)
                obj.p1(pos) = obj.p1in - sum(deltap_i(1:pos));  
            end
            obj.delta_p1 = deltap_i;

        end
        %% HEAT TRANSFER COEFFICIENT OUTSIDE
        function alpha_out(obj)
            persistent lambda40 baseAlpha
            if isempty(lambda40)
                lambda40=tppAir('lambda',40+273.15)^0.6;

                d_sv=120e-6;
                baseAlpha=1200.*d_sv.^-0.36./150e-6.^-0.36;
            end
            alpha=baseAlpha.*tppAir('lambda',obj.T3).^0.6./lambda40;
            
            obj.alpha_a = sum(alpha)/obj.nCells;
        end    

        %% HEAT EXCHANGER CALCULATION
        % DESCRIPTIVE TEXT 
        function hex_calculation(obj,idx)
            da = obj.di + 2*obj.e;
            A1=repmat(da.*pi.*obj.N*obj.l/obj.nCells,1,obj.nCells);         
            
            if obj.charge
                T1in=obj.T1left;
                T3in=obj.T3right;
                h1in=obj.h1left;
                h3in=obj.h3right;
            else
                T1in=fliplr(obj.T1right);
                T3in=fliplr(obj.T3left);
                h1in=fliplr(obj.h1right);
                h3in=fliplr(obj.h3left);
                
                obj.T1=fliplr(obj.T1);
                obj.T3=fliplr(obj.T3);
                
                obj.p1=fliplr(obj.p1);
                obj.rho1=fliplr(obj.rho1);
            end
          
            err=repmat(10,1,5);
            QDot=0;
            count=0;
            isConv=false;
            while err(end)>1 || ~isConv
                
                obj.alpha_in();
                obj.alpha_out();
                obj.k = ((1./obj.alpha_a)+ da./(2.*obj.lambda_st).*log(da./obj.di)+da./(obj.alpha_i.*obj.di)).^-1;

                dT = obj.T1 -obj.T3;
                QDot1=obj.F.*obj.k.*A1.*dT; %Charge / Discharge
                QDot3=QDot1;
                
                %energy balance                
                
                h1out=h1in-QDot1./obj.mDot1;
                h1in(2:end)=h1out(1:end-1);
                
                T1out=IF97.T_ph(obj.p1,h1out);
                
                %Consistency check
                if obj.charge
                    check=T1out<T3in;
                else
                    check=T1out>T3in;    
                end
                T1out(check)=T3in(check);
                
                T1in(2:end)=T1out(1:end-1);
                obj.T1=mean([T1in;T1out],1);
               
                h3out=h3in+QDot3./obj.mDot3;
                h3in(1:end-1)=h3out(2:end);
                
                T3out=SiO2.T_h(h3out);
                
                %Consistency check
                if obj.charge
                    check=T3out>T1in;
                else
                    check=T3out<T1in;    
                end
                T3out(check)=T3out(check);
                
                T3in(1:end-1)=T3out(2:end);
                obj.T3=mean([T3in;T3out],1);
                
                obj.deltap_i();
                obj.rho1=hex.rho(obj.p1,mean([h1in;h1out],1));
                
                %Abbruchbedingung
                QDotsum=sum(QDot3);
                err=circshift(err,-1);
                err(end)=abs(QDotsum-QDot);
                QDot=QDotsum;
                obj.Qdot = QDot;
           

                %Ausgabe-Information
                if mod(count,10)==0
                    switch idx
                        case 1
                            disp(['count=',num2str(count),', err=',num2str(err(end)),' QDot=',num2str(QDot)]);
                        case 2
                            disp(['count=',num2str(count),', err=',num2str(err(end)),' TH2Oout=',num2str(obj.T1(end))]);    
                    end
                end

                %Notfall-Abbruch ab 10000 Iterationen
                count=count+1;
                if count>10000
                    break
                end

                %Konvergenz-Bedingung
                derr=diff(err);
                isConv=(all(derr<0) && err(end)/err(1)>0.7) || max(abs(derr))<0.2;

            end
                     
            if obj.charge
                obj.T1left=T1in;
                obj.T1right=T1out;
                
                obj.T3left=T3out;
                obj.T3right=T3in;
                
                obj.h1left=h1in;
                obj.h1right=h1out;
                
                obj.h3left=h3out;
                obj.h3right=h3in;
            else
                obj.T1left=fliplr(T1out);
                obj.T1right=fliplr(T1in);
                
                obj.T3left=fliplr(T3in);
                obj.T3right=fliplr(T3out);
                
                obj.h1left=fliplr(h1out);
                obj.h1right=fliplr(h1in);
                
                obj.h3left=fliplr(h3in);
                obj.h3right=fliplr(h3out);
                
                obj.T1=fliplr(obj.T1);
                obj.T3=fliplr(obj.T3);
                
                obj.p1=fliplr(obj.p1);
                obj.delta_p1=fliplr(obj.delta_p1);
                obj.rho1=fliplr(obj.rho1);
            end
        obj.count_a = count;
        obj.err_a = err;
        end
        
        %% Tx DIAGRAM
        function Tx_diagram(obj)
            x=linspace(obj.l./obj.nCells,obj.l,obj.nCells);
            x(x<obj.di) = obj.di;
            plot(x,[obj.T1;obj.T3]);
            legend('Wasser','Sand')
            xlabel("x in m")
            ylabel("T in K")
        end
        %% TQ DIAGRAM NEW
        % plots a new TQ diagram 
        function plot_TQ(obj)
            Qdot_n=linspace(0,abs(obj.Qdot),obj.nCells+1);  
            if obj.T1(1) > obj.T3(1) % charge
                plot(Qdot_n,[[obj.T1left,obj.T1right(end)];[obj.T3left(1),obj.T3right]]);
            else % discharge
                plot(Qdot_n,[[obj.T1left(1),obj.T1right];[obj.T3left,obj.T3right(end)]]);
            end
            legend('Wasser','Sand')
            xlabel("QDot in W")
            ylabel("T in K")
        end

        %% SECTION TITLE
        % px DIAGRAM
        function px_diagram(obj)
            x=linspace(obj.l./obj.nCells,obj.l,obj.nCells);
            x(x<obj.di) = obj.di;

            plot(x,obj.p1);
            xlabel("x in m")
            ylabel("p in Pa")
           
        end 
    end
    
    % density characteristic
    methods (Static)
        function rho=rho(p,h)
            hsatl=IF97.h(p,NaN,0);
            hsatv=IF97.h(p,NaN,1);
            is2phase=hsatl<=h & h<=hsatv;
            
            x=NaN(size(h));
            x(is2phase)=(h(is2phase)-hsatl(is2phase))./(hsatv(is2phase)-hsatl(is2phase));
            
            T=NaN(size(h));
            T(~is2phase)=IF97.T_ph(p(~is2phase),h(~is2phase));
            
            rho=IF97.v(p,T,x).^-1;
        end
    end
 end




