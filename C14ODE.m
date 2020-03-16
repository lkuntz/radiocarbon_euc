function [dCdt,y] = C14ODE(t,x,u,ft,f,vt,v,vot,vo,At,A,C_ab,C_abt,C_euc,C_euct,frac_a,frac_d,C_d,m,n,a,ha,ho,voff1,voff2,voff3,varargin)
    y = x;
    frac_euc = 1-frac_d-frac_a;
    h = ha*(-sin(pi*[1:length(ft)]/6))+ho;
    f = f./h; f = f(:);
    f = a*interp1(ft+voff3,f,t);
    A = interp1(At,A,t);
    v = interp1(vt+voff1,v,t);
    vo = interp1(vot+voff2,vo,t);
    C_ab = interp1(C_abt,C_ab,t);
    C_euc = interp1(C_euct,C_euc,t);
    V = m*v*heaviside(t-18*12)+n*v*heaviside(18*12+eps-t);
    Vo = m*vo*heaviside(t-18*12)+n*vo*heaviside(18*12+eps-t);
    
    dCdt = -x*(1.3*f+Vo)+f*A+V*(frac_euc*C_euc+frac_d*C_d+frac_a*C_ab);

% (t,C,u,ft,f,vt,v,C_ab,C_abt,C_euc,C_euct,frac_a,frac_d,m,n,C_d,varargin)
%     frac_euc = 1-frac_d-frac_a;
%     A = interp1(At,A,t);
%     v = interp1(vt,v,t);
%     f = interp1(ft,f,t);
%     C_ab = interp1(C_abt,C_ab,t);
%     C_euc = interp1(C_euct,C_euc,t);
%     V = m*v*heaviside(t-18*12)+n*v*heaviside(18*12+eps-t);
%     
%     dCdt = -C*(1.3*f+v)+f*A+v*(frac_euc*C_euc+frac_d*C_d+frac_a*C_ab);
end