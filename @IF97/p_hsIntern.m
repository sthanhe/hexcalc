function p=p_hsIntern(h,s,regions,region2,region3)
    persistent p_sat273
    if isempty(p_sat273)
        p_sat273=IF97.p_satIntern(273.15);
    end
    
    
    p=NaN(size(h));
    p(regions(1))=p_hs1(h(regions(1)),s(regions(1)));
    
    p(region2(1))=p_hs2a(h(region2(1)),s(region2(1)));
    p(region2(2))=p_hs2b(h(region2(2)),s(region2(2)));
    p(region2(3))=p_hs2c(h(region2(3)),s(region2(3)));
    
    p(region3(1))=p_hs3a(h(region3(1)),s(region3(1)));
    p(region3(2))=p_hs3b(h(region3(2)),s(region3(2)));
    
%     for i=find(regions(4))
%         p(i)=fzero(@(p) deltax(h(i),s(i),p),[p_sat273,IF97.p_c-10^-10]);
%     end
end


function p=p_hs1(h,s)
    persistent I J n p_star h_star s_star
    if isempty(h_star)
        vars=load('Region1constants.mat','I_hs','J_hs','n_hs');
        I=vars.I_hs;
        J=vars.J_hs;
        n=vars.n_hs;
        p_star=100*10^6;
        h_star=3400*10^3;
        s_star=7.6*10^3;
    end
    
    
    if ~isempty(h)
        p=p_star.*sum(n.*(h./h_star+0.05).^I.*(s./s_star+0.05).^J,1);
    else
        p=NaN(1,0);
    end
end


function p=p_hs2a(h,s)
    persistent I J n p_star h_star s_star
    if isempty(h_star)
        vars=load('Region2constants.mat','I_hsa','J_hsa','n_hsa');
        I=vars.I_hsa;
        J=vars.J_hsa;
        n=vars.n_hsa;
        p_star=4*10^6;
        h_star=4200*10^3;
        s_star=12*10^3;
    end
    
    
    if ~isempty(h)
        p=p_star.*sum(n.*(h./h_star-0.5).^I.*(s./s_star-1.2).^J,1).^4;
    else
        p=NaN(1,0);
    end
end


function p=p_hs2b(h,s)
    persistent I J n p_star h_star s_star
    if isempty(h_star)
        vars=load('Region2constants.mat','I_hsb','J_hsb','n_hsb');
        I=vars.I_hsb;
        J=vars.J_hsb;
        n=vars.n_hsb;
        p_star=100*10^6;
        h_star=4100*10^3;
        s_star=7.9*10^3;
    end
    
    
    if ~isempty(h)
        p=p_star.*sum(n.*(h./h_star-0.6).^I.*(s./s_star-1.01).^J,1).^4;
    else
        p=NaN(1,0);
    end
end


function p=p_hs2c(h,s)
    persistent I J n p_star h_star s_star
    if isempty(h_star)
        vars=load('Region2constants.mat','I_hsc','J_hsc','n_hsc');
        I=vars.I_hsc;
        J=vars.J_hsc;
        n=vars.n_hsc;
        p_star=100*10^6;
        h_star=3500*10^3;
        s_star=5.9*10^3;
    end
    
    
    if ~isempty(h)
        p=p_star.*sum(n.*(h./h_star-0.7).^I.*(s./s_star-1.1).^J,1).^4;
    else
        p=NaN(1,0);
    end
end


function p=p_hs3a(h,s)
    persistent I J n p_star h_star s_star
    if isempty(h_star)
        vars=load('Region3constants.mat','I_hsa','J_hsa','n_hsa');
        I=vars.I_hsa;
        J=vars.J_hsa;
        n=vars.n_hsa;
        p_star=99*10^6;
        h_star=2300*10^3;
        s_star=4.4*10^3;
    end
    
    
    if ~isempty(h)
        p=p_star.*sum(n.*(h./h_star-1.01).^I.*(s./s_star-0.75).^J,1);
    else
        p=NaN(1,0);
    end
end


function p=p_hs3b(h,s)
    persistent I J n p_star h_star s_star
    if isempty(h_star)
        vars=load('Region3constants.mat','I_hsb','J_hsb','n_hsb');
        I=vars.I_hsb;
        J=vars.J_hsb;
        n=vars.n_hsb;
        p_star=16.6*10^6;
        h_star=2800*10^3;
        s_star=5.3*10^3;
    end
    
    
    if ~isempty(h)
        p=p_star./sum(n.*(h./h_star-0.681).^I.*(s./s_star-0.792).^J,1);
    else
        p=NaN(1,0);
    end
end


% function deltax=deltax(h,s,p)
%     T=IF97.T_satIntern(p);
%     
%     if T<=623.15
%         h0=IF97.h_pTx(p,T,0,@(x) x==1,true);
%         h1=IF97.h_pTx(p,T,1,@(x) x==2,true);
% 
%         s0=IF97.s_pTx(p,T,0,@(x) x==1,true);
%         s1=IF97.s_pTx(p,T,1,@(x) x==2,true);
%     else
%         rho0=IF97.rho_satLintern(T);
%         rho1=IF97.rho_satVintern(T);
%         
%         h0=IF97.h_rhoTintern(rho0,T,@(x) x==3,rho0,rho1,0);
%         h1=IF97.h_rhoTintern(rho1,T,@(x) x==3,rho0,rho1,1);
% 
%         s0=IF97.s_rhoTintern(rho0,T,@(x) x==3,rho0,rho1,0);
%         s1=IF97.s_rhoTintern(rho1,T,@(x) x==3,rho0,rho1,1);
%     end
% 
%     deltax=(h-h0)/(h1-h0)-(s-s0)/(s1-s0);
% end




