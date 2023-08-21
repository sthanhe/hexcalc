classdef recuperator < handle

    % VDI Wärmeatlas - C1 Berechnung von Wärmeübertragern S.37
    % Reiner Kreuzstrom

    %1: exhaust air - hot side
    %1: fluidization air - could side

    properties
        %input
        Th_in;
        Tc_in;
        mDot4;
        eps %Epsilon (P) dimensionless temperature change

        %output
        Tc_out;
        Th_out;
        QDot;
        kA;
        NTU1;
        NTU2;
    end

    methods
        function obj = recuperator(Th_in, Tc_in, mDot4)
            if nargin==0
                return
            end
            obj.Th_in = Th_in;
            obj.Tc_in = Tc_in;
            obj.mDot4 = mDot4;          
        end

        function desing(obj,p,eps)
            if nargin>2
                obj.eps = eps;
            end
          
            % obj.Tc_out = obj.eps*(obj.Th_in-obj.Tc_in)+obj.Tc_in; %Formel (10) VDI
            obj.Tc_out = obj.Th_in - 20; %Grädigkeit
            % p = 1.8*10^5;

            hh_in = DryAir.h(p,obj.Th_in);
            hc_in = DryAir.h(p,obj.Tc_in);
            hc_out = DryAir.h(p,obj.Tc_out);

           
            obj.QDot = -obj.mDot4*(hc_in-hc_out);
            hh_out = hh_in - obj.QDot/obj.mDot4;
            obj.Th_out = DryAir.T_ph(p,hh_out); 

            dTgr=obj.Th_in-obj.Tc_out;
            dTkl=obj.Th_out-obj.Tc_in;

            LMTD = (dTgr-dTkl)./log(dTgr/dTkl);
%             Theta = LMTD/(Th_in - Tc_in);

            NTU2g = (obj.Tc_out-obj.Tc_in)/LMTD;
            NTU1g = (obj.Th_in-obj.Th_out)/LMTD;

            % Wärmekapazitätsstrom
            WDot1=obj.mDot4*(hh_in-hh_out)/(obj.Th_in-obj.Th_out);
            WDot2=obj.mDot4*(hc_in-hc_out)/(obj.Tc_in-obj.Tc_out);

            %umrechnung in NTU REAL
            R1 = WDot1/WDot2;
%             R2 = 1/R1;

            a = 0.433; 
            b = 0.16;
            c = 0.267;
            d = 0.5;

            dT_NTU = 1;
            counter = 0;
            obj.NTU1 = NTU1g;
            while dT_NTU > 1e-6 && counter<1000
                NTU1_old = obj.NTU1;
                obj.NTU1 = NTU1g*(1+a*R1^(d*b)*NTU1_old^b)^c;
                dT_NTU = abs(NTU1_old-obj.NTU1);
                counter = counter + 1;
            end

            F = NTU1g/obj.NTU1;
            obj.NTU2 = NTU2g/F;

            obj.kA = obj.NTU1*WDot1;

        end

        function recalculate(rec_r,rec_d,mDot4,Th_in,Tc_in,p)
            % Pi=Pi(NTUi,Ri)
            % p = 1.8*10^5;
            rec_r.NTU1 = rec_d.NTU1;
            rec_r.NTU2 = rec_d.NTU2;
            rec_r.mDot4 = mDot4;

            if nargin < 4
                rec_r.Th_in = rec_d.Th_in;
                rec_r.Tc_in = rec_d.Tc_in;
            else
                rec_r.Th_in = Th_in;
                rec_r.Tc_in = Tc_in;
            end

            hh_in = DryAir.h(p,rec_r.Th_in);
%             hc_in = DryAir.h(rec_r.Tc_in);

            %assumption
            % rec_r.Th_out = rec_r.Th_in-(rec_d.Th_in-rec_d.Th_out);
            % rec_r.Tc_out = rec_r.Tc_in+(rec_d.Tc_out-rec_d.Tc_in);
            rec_r.Tc_out = rec_r.Th_in - 20 - 273.15; %Grädigkeit 25°C
            hc_in = DryAir.h(p,rec_r.Tc_in);
            hc_out = DryAir.h(p,rec_r.Tc_out);
            rec_r.QDot = abs(rec_r.mDot4*(hc_out-hc_in));
            hh_out = hh_in - rec_r.QDot/rec_r.mDot4;
            rec_r.Th_out = DryAir.T_ph(p,hh_out); 


            % hh_out = DryAir.h(p,rec_r.Th_out);
%             hc_out = DryAir.h(rec_r.Tc_out); 

            dT=1;
            dTc=1;
            counter=0;
            % aktuell wir keine while Schleife benötigt
            while dT>1e-6 && dTc>1e-6 && counter<1000 
                Th_old = rec_r.Th_out;
                Tc_old = rec_r.Tc_out;

                R2 = rec_r.NTU1/rec_r.NTU2;
                R1 = 1/R2;

                %correction factor
                a = 0.433; % kann noch mit Tabelle und automatischer Auswahl verbessert werden.
                b = 0.16;
                c = 0.267;
                d = 0.5;

                F = 1/(1+a*R1^(d*b)*rec_r.NTU1^b)^c;

                P1 = (1 - exp((R1-1)*rec_r.NTU1*F))/(1-R1*exp((R1-1)*rec_r.NTU1*F));
                P2 = (1 - exp((R2-1)*rec_r.NTU2*F))/(1-R2*exp((R2-1)*rec_r.NTU2*F));              
                
                rec_r.Th_out = rec_r.Th_in - P1*(rec_r.Th_in - rec_r.Tc_in);
                rec_r.Tc_out = rec_r.Tc_in + P2*(rec_r.Th_in - rec_r.Tc_in);

                dT = abs(Th_old-rec_r.Th_out);
                dTc = abs(Tc_old-rec_r.Tc_out);
                counter = counter + 1;
            end
          
            WDot1=rec_r.mDot4*(hh_in-hh_out)/(rec_r.Th_in-rec_r.Th_out);
%             WDot2=rec_r.mDot4*(hc_in-hc_out)/(rec_r.Tc_in-rec_r.Tc_out); 
            rec_r.kA = rec_r.NTU1 * WDot1;
           
            rec_r.QDot = WDot1*(rec_r.Th_in-rec_r.Th_out);
        end
    end
end

                




                
                

            
    
  
