function [ d14C ] = Delta14Todelta( Delta14C,d13C )
%Delta14ToConcentration takes in the measured Delta14C value and Total DIC
%concentration and returns the concentration of 14C isotope
if ~exist('d13C','var')
    d13C = -19;
end
Rstd = 1/8267; %absolute international standard activity

%d14C = (Delta14C+2/1000*(d13C+25))/(1-2/1000*(d13C+25));
d14C = (Delta14C+2*(d13C+25))/(1-2/1000*(d13C+25));


end

