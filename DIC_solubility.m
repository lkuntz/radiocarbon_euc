function [Ko,k] = DIC_solubility(Temp_C,Salt_ppt,wind)
%DIC_SOLUBILITY determines the gas transfer velocity (k) and the solubility
%of C (Ko) based on the temperature, salinity, and wind

Temp_K = Temp_C+273.15;

%Coefficients from Weiss, 1974
A1 = -58.0931;
A2 = 90.5069;
A3 = 22.2950;
B1 = 0.027766;
B2 = -.025888;
B3 = 0.0050578;

%in moles/L/atm
Ko = exp(A1+A2*(100/Temp_K)+A3*log(Temp_K/100)...
    +Salt_ppt*(B1+B2*(Temp_K/100)+B3*(Temp_K/100)^2));

Ko = Ko*1000; %in moles/m^3/atm

%From Wang et al, 2006
Sc = 2073.1-125.62*Temp_C+3.6276*Temp_C^2-0.043219*Temp_C^3;
k = 0.251*wind^2*(Sc/660)^(-1/2); %in cm/hr
k = k*24*365.25/12/100; %in m/month
end

