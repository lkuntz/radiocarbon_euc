function [ Delta14C ] = ConcentrationToDelta14( C14,DIC )
%ConcentrationToDelta14 takes in the the concentraion of 14C isotopes and
%total DIC concentration and returns the Delta14C
d13C = -19;
Rstd = 1/8267; %absolute international standard activity

d14C = ((C14./DIC)/Rstd-1)*1000;

Delta14C = d14C-2.*(d13C+25).*(1+d14C/1000);

end

