classdef compressor < handle 
    properties
        % compressor
        pin = 10^5                % ambient pressure (Pa)
        pout;                     % pressure after compression (Pa)
        dp_c;                     % pressure difference compressor (Pa)
        Tin = 25 + 273.15         % temperature ambient air (K)
        Tout;                     % temperature after compression (K)
        eta = 0.7                 % efficiency
        wt;                       % technical work of the compressor (J/kg)

        % sand bed
        h_bed;                    % height of the sand bed (m)
        A_bed;                    % area of the bed (m^2)
        delta_p_bed;              % pressure loss of the bed (Pa)
        p_bed = 0;                % pressure of the bed on top (Pa)
        dp = 120 * 10^-6;         % diameter of the sand particle (m)
        rho_p = 2650;             % density of the particle (kg/m^3)
        eps_por = 0.5;            % porosity of the sand bed (-)
    end

    methods

        function obj = compressor()
            if nargin==0
                return
            else

            end
        end

        function calculate(obj,factor)
            obj.delta_p_bed = FluBed.deltaP(obj.h_bed,obj.eps_por,obj.rho_p);
            obj.dp_c = compressor.deltap(obj.delta_p_bed ,factor);
            obj.pout = obj.pin + obj.dp_c;
            obj.Tout = compressor.Tout_dp(obj.Tin, obj.eta, obj.pin, obj.pout);
            obj.wt = compressor.w_t(obj.pin, obj.pout, obj.Tin, obj.eta);
            obj.p_bed = compressor.pbed(obj.pout, obj.delta_p_bed);
        end
    end

    methods(Static)
        function T2 = Tout_dp(T1, eta, p1, p2)
            T2 = T1 + eta.^-1.*T1.*((p2/p1).^(0.4/1.4)-1);
        end

        function dp_c = deltap(dp_bed, factor)
            dp_c = dp_bed * factor;
        end

        function w = w_t(p1, p2, T1, eta)
            cp = DryAir.c_p(p1,T1);
            w = eta^-1*cp*T1*((p2/p1)^(0.4/1.4)-1);
        end

        function p = pbed(pout, dp_bed)
            p = pout - dp_bed;
        end

    end
end