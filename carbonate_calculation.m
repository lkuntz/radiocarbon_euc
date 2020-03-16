function [pCO2] = carbonate_calculation(Temp_C,DIC)
%Calculates the partial pressure of CO2 using Temperature and DIC as inputs

Temp_K = Temp_C+273.15;
DIC_molperkg = DIC/1025; %[mol/kg]

S = 34.78;
Boron = 1.179e-5*S;

Alk = 65.89*S+4.124; %[microunits/kg^3]
Alk = Alk*10^-6; %[units/kg]

K0 = exp(-60.2409 + 9345.17/Temp_K + 23.3585*log(Temp_K/100) ...
          + S * (0.023517 - 0.00023656*Temp_K +0.0047036*(Temp_K/100)^2) );
      
K1 = exp(2.18867 - 2275.036/Temp_K - 1.468591 * log(Temp_K) ...
          + (-0.138681 - 9.33291/Temp_K) * sqrt(S) + 0.0726483*S...    
          - 0.00574938 * S ^1.5);

K2 = exp(-0.84226 - 3741.1288/Temp_K -1.437139 * log(Temp_K)...
          + (-0.128417 - 24.41239/Temp_K)*sqrt(S) + 0.1195308 * S ...  
          - 0.0091284 * S ^1.5 );

Kb = exp( (-8966.90 - 2890.51*sqrt(S) - 77.942*S ...
            + 1.726 * S ^1.5 - 0.0993*S^2) / Temp_K ...               
            + (148.0248 + 137.194 * sqrt(S) + 1.62247 * S) ...      
            + (-24.4344 - 25.085 * sqrt(S) - 0.2474 * S) * log(Temp_K)...
            + 0.053105 * sqrt(S) * Temp_K);
        
pH = 8;
H = 10e-8;
diff_H = H;

iter = 1;
while diff_H>eps
    H_old = H;                      %remember old value of H

  %solve Tans' equation 13 for carbonate alkalinity from TA
  CA = Alk - (Kb/(Kb+H)) * Boron;      

  % solve quadratic for H (Tans' equation 12)
  a = CA;
  b = K1 * (CA - DIC_molperkg);
  c = K1 * K2 * (CA - 2 * DIC_molperkg);
  H = (-b + sqrt(b^2 - 4. * a * c) ) / (2. * a);  

  % How different is new estimate from previous one?
  diff_H = abs(H - H_old);
  iter = iter + 1;
end

CO2aq = CA / (K1/H + 2*K1*K2/H^2);
pCO2 = CO2aq / K0 * 1.e6 ;

end


