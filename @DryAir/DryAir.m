%% Property functions of dry air
%GNU General Public License v3.0
%CC-By Stefan Thanheiser: https://orcid.org/0000-0003-2765-1156
%
%
%This class describes the thermo-physical properties of dry air based on a
%lookup-table according to:
%
%Span, R. (2019). D2.2 Thermophysikalische Stoffwerte von trockener Luft. 
%In: Stephan, P., Kabelac, S., Kind, M., Mewes, D., Schaber, K., Wetzel, 
%T. (eds) VDI-Wärmeatlas. Springer Reference Technik(). Springer Vieweg, 
%Berlin, Heidelberg. https://doi.org/10.1007/978-3-662-52989-8_13
%
%
%Requires all files packaged in the class folder and on the MATLAB path
%
%Required products:
%   - MATLAB, version 9.14
%   - Curve Fitting Toolbox, version 3.9
%Data files:
%   - dryAirFits.mat
%Additional classes:
%   - implExp


classdef DryAir
    %All parameters and results in SI base units
    
    %% Constants
    properties(Constant)
        M=28.9583e-3;   %molar mass
        R=287.12;       %specific gas constant
    end
    
    
    %% prop(p,T) functions
    % folgende Funktionen werden nur als Platzhalter mit NaN ausgegeben
    methods(Static)
        function i2ph = is2phase(Rho,T)
            i2ph =  Rho == T;
        end

        function rho_sat = rho_satL(T)
            rho_sat = NaN(length(T));
        end

        function rho_sat = rho_satV(T)
            rho_sat = NaN(length(T));
        end

        function xDot = x(~,T)
            xDot = NaN(length(T));
        end

        function s = sigma(T)
            s = NaN(length(T));
        end
    end
    methods(Static)
        function rho=rho(p,T)
            %Density
            persistent fx
            if isempty(fx)
                fxStruct=load('@DryAir\dryAirFits.mat','rho');
                fx=fxStruct.rho;
            end
            
            
            rho=fx(p,T);
        end
        
        
        function h=h(p,T)
            %Specific enthalpy
            persistent fx
            if isempty(fx)
                fxStruct=load('@DryAir\dryAirFits.mat','h');
                fx=fxStruct.h;
            end
            
            h=fx(p,T);
        end
        
        
        function s=s(p,T)
            %Specific entropy
            persistent fx
            if isempty(fx)
                fxStruct=load('@DryAir\dryAirFits.mat','s');
                fx=fxStruct.s;
            end
            
            s=fx(p,T);
        end
        
        
        function c_p=c_p(p,T)
            %Specific isobaric heat capacity
            persistent fx
            if isempty(fx)
                fxStruct=load('@DryAir\dryAirFits.mat','c_p');
                fx=fxStruct.c_p;
            end
            
            c_p=fx(p,T);
        end
        
        
        function c_v=c_v(p,T)
            %Specific isochoric heat capacity
            persistent fx
            if isempty(fx)
                fxStruct=load('@DryAir\dryAirFits.mat','c_v');
                fx=fxStruct.c_v;
            end
            
            c_v=fx(p,T);
        end
        
        
        function kappa=kappa(p,T)
            %Isentropic exponent
            kappa=DryAir.c_p(p,T)./DryAir.c_v(p,T);
        end
        
        
        function beta=beta(p,T)
            %Coefficient of thermal expansion
            persistent fx
            if isempty(fx)
                fxStruct=load('@DryAir\dryAirFits.mat','beta');
                fx=fxStruct.beta;
            end
            
            beta=fx(p,T);
        end
        
        
        function w_s=w_s(p,T)
            %Speed of sound
            persistent fx
            if isempty(fx)
                fxStruct=load('@DryAir\dryAirFits.mat','w');
                fx=fxStruct.w;
            end
            
            w_s=fx(p,T);
        end
        
        
        function lambda=lambda(p,T)
            %Thermal conductivity
            persistent fx
            if isempty(fx)
                fxStruct=load('@DryAir\dryAirFits.mat','lambda');
                fx=fxStruct.lambda;
            end
            
            lambda=fx(p,T);
        end
        
        
        function eta=my(p,T)
            %Dynamic viscosity
            persistent fx
            if isempty(fx)
                fxStruct=load('@DryAir\dryAirFits.mat','eta');
                fx=fxStruct.eta;
            end
            
            eta=fx(p,T);
        end
        
        
        function ny=ny(p,T)
            %Kinematic viscosity
            persistent fx
            if isempty(fx)
                fxStruct=load('@DryAir\dryAirFits.mat','ny');
                fx=fxStruct.ny;
            end
            
            ny=fx(p,T);
        end
        
        
        function a=a(p,T)
            %Thermal diffusivity            
            persistent fx
            if isempty(fx)
                fxStruct=load('@DryAir\dryAirFits.mat','a');
                fx=fxStruct.a;
            end
            
            a=fx(p,T);
        end
        
        
        function Pr=Pr(p,T)
            %Prandtl number
            persistent fx
            if isempty(fx)
                fxStruct=load('@DryAir\dryAirFits.mat','Pr');
                fx=fxStruct.Pr;
            end
            
            Pr=fx(p,T);
        end

        % @TB
        function v = v(p,T)
            v = DryAir.rho(p,T).^-1;
        end
    end
    
    
    %% Other property functions
    methods(Static)
        function T=T_ph(p,h)
            %Backwards-equation for temperature as function of pressure and 
            %specific enthalpy
            persistent hbounds
            if isempty(hbounds)
                hbounds=[-331200,1152000];
            end

            sz=implExp.size(p,h);
            [p,h]=implExp.normalize(sz,p,h);

            T=NaN(sz);
            for i=1:numel(T)
                if h(i)>=hbounds(1) && h(i)<=hbounds(2)
                    T(i)=fzero(@(T) DryAir.h(p(i),T)-h(i),hbounds);
                end
            end
        end

        % @TB
        function p = p(rho,T)
            p = NaN(1,size(T,2));
            for i=1:size(T,2) 
                p(i) = fzero(@(p) DryAir.rho(p,T(i))-rho(i),[0, 10^7]);
            end
        end

        function h = h_rhoT(rho,T)
            p = DryAir.p(rho,T);
            h = DryAir.h(p,T);
        end 

    end

    % @TB
    % Werte für Stickstoff
    properties(Constant)
        cf = 0.8;    %properties of the fluid                                                    
        q0 = 10000;  %substance-specific values
        d0 = 10^-2 ;
        Ra0 = 10^-6;
        alpha_0 = 4380; %W/m^2K
        p_c = 34 * 10^5 %Pa
    end
    
    
    %% Auxilliary Functions
    methods(Static)
        function createConstants()
            %Creates the lookup-constants from the Excel-file
            clear('DryAir');


            tab=readtable('dryAir.xls','Sheet','h');
            T=tab.Properties.VariableNames(2:end);
            T=strrep(T,'T','');
            T=str2double(T)+0.15;
            p=tab.p.*10^5;

            
            tab=readtable('dryAir.xls','Sheet','h');
            vals=tab{:,2:end}'.*10^3;
            h=DryAir.fit_h(p,T,vals);

            tab=readtable('dryAir.xls','Sheet','rho');
            vals=tab{:,2:end}';
            rho=DryAir.fit_rho(p,T,vals);

            tab=readtable('dryAir.xls','Sheet','s');
            vals=tab{:,2:end}'.*10^3;
            s=DryAir.fit_s(p,T,vals);

            tab=readtable('dryAir.xls','Sheet','c_p');
            vals=tab{:,2:end}'.*10^3;
            c_p=DryAir.fit_c_p(p,T,vals);

            tab=readtable('dryAir.xls','Sheet','c_v');
            vals=tab{:,2:end}'.*10^3;
            c_v=DryAir.fit_c_v(p,T,vals);

            tab=readtable('dryAir.xls','Sheet','beta');
            vals=tab{:,2:end}'.*10^-3;
            beta=DryAir.fit_beta(p,T,vals);

            tab=readtable('dryAir.xls','Sheet','w');
            vals=tab{:,2:end}';
            w=DryAir.fit_w(p,T,vals);

            tab=readtable('dryAir.xls','Sheet','lambda');
            vals=tab{:,2:end}'.*10^-3;
            lambda=DryAir.fit_lambda(p,T,vals);

            tab=readtable('dryAir.xls','Sheet','eta');
            vals=tab{:,2:end}'.*10^-6;
            eta=DryAir.fit_eta(p,T,vals);

            tab=readtable('dryAir.xls','Sheet','ny');
            vals=tab{:,2:end}'.*10^-7;
            ny=DryAir.fit_ny(p,T,vals);

            tab=readtable('dryAir.xls','Sheet','a');
            vals=tab{:,2:end}'.*10^-7;
            a=DryAir.fit_a(p,T,vals);

            tab=readtable('dryAir.xls','Sheet','Pr');
            vals=tab{:,2:end}';
            Pr=DryAir.fit_Pr(p,T,vals);


            save('@DryAir\dryAirFits.mat','h','rho','s','c_p','c_v','beta','w','lambda','eta','ny','a','Pr');
        end
    end


    methods(Static,Access=private)
        [fitresult, gof] = fit_h(p, T, h)
        [fitresult, gof] = fit_rho(p, T, rho)
        [fitresult, gof] = fit_w(p, T, w)
        [fitresult, gof] = fit_s(p, T, s)
        [fitresult, gof] = fit_ny(p, T, ny)
        [fitresult, gof] = fit_Pr(p, T, Pr)
        [fitresult, gof] = fit_lambda(p, T, lambda)
        [fitresult, gof] = fit_eta(p, T, eta)
        [fitresult, gof] = fit_c_v(p, T, c_v)
        [fitresult, gof] = fit_c_p(p, T, c_p)
        [fitresult, gof] = fit_beta(p, T, beta)
        [fitresult, gof] = fit_a(p, T, a)
    end
end




