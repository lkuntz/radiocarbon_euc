function [ C14 ] = Delta14ToConcentration( Delta14C,DIC )
%Delta14ToConcentration takes in the measured Delta14C value and Total DIC
%concentration and returns the concentration of 14C isotope
d13C = -19;
Rstd = 1/8267; %absolute international standard activity

%d14C = (Delta14C+2/1000*(d13C+25))/(1-2/1000*(d13C+25));
d14C = (Delta14C+2*(d13C+25))/(1-2/1000*(d13C+25));

C14 = DIC*Rstd*(d14C/1000+1);

end

