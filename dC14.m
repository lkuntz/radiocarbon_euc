function [ Cp ] = dC14( a,x )
%Calculate the time-series of C14 concentration predicted using the
%measured values in the N and S pacific as well as the wind data (proxy for 
%upwelling magnitude variations), atmospheric data (penetration),
%horizontal advection of northern waters, and sea level in the west pacific
%(proxy for amount of deep water upwelling)

Cn = x(:,1); %northern C14 record
Cs = x(:,2); %southern C14
Cd = Delta14ToConcentration(-100,1); %assume deep water of -100 per mil C14
Catm = ConcentrationToDelta14(x(:,3),1); %atmospheric data as a delta
Cv = x(:,4); %off equatorial east pacific C14
Cp = x(:,5); %record of C14 in east Pacific (only actually using the first value)
Vup = a(4)*sin([0:307]*pi/6)+a(3); %upwelling strength varies seasonally

%Currently neither of these are used, but have been used in older versions
W = x(:,6); %wind data (NOT USED)
S = x(:,7); %western sea level data (NOT USED)

%different possibilities for meridional advection parameterization
%V = a(5)*cos([0:307]*pi/6)+a(6); %seasonal cycle
V = a(5); %constant
%V=V(:);

for i=1:length(Cn)


    DCp = ConcentrationToDelta14(Cp(i),1);
    
    %Flux into east pacific box is a mixture of atmospheric penetration
    %(proportional to the difference in atmosphere and previous C14
    %concentration in east Pacific), upwelled deep water, upwelled north
    %and south pacific water, and meridionally advected water.  By mass
    %continuity, the amount of water added equals the amount removed
    %(amount removed is taken to have C14 of the previous time step)
    F(i) = a(1)*(Catm(i)-DCp)+a(2)*Vup(i)*(Cd-Cp(i)) + Vup(i)*(1/3*Cn(i)+2/3*Cs(i)-Cp(i))...
        +V*(Cv(i)-Cp(i));
    
    %Concentration of C14 in the next time step is taken by adding the flux
    Cp(i+1) = Cp(i)+F(i);
end

F=double(F);

Cp = double(Cp)*1e4;

end


%two previous used version of the flux parameterization - the difference
%stems from the fact that the wind data and sea level pressure data are
%added to try and account for variations in upwelling magnitude and
%N/S/deep water mixture of waters.  The first equation only includes wind
%and sea level pressure as proxies for amount of upwelling and ratio, while
%the second equation uses them to describe variations from the seasonal
%climatological cycle
% F = (1/3*Cn.*(a(1)*W-a(2)*S)+2/3*Cs.*(a(1)*W-a(2)*S)+Cd*a(2)*S+V.*Cv-...
%     Cp.*(a(1)*W+a(2)*S)+a(3)*Catm);

%     F(i) = (1/3*Cn(i)*(a(1)*W(i)+Vup(i)-a(2)*S(i))+2/3*Cs(i).*(a(1)*W(i)+Vup(i)-a(2)*S(i))+Cd*a(2)*S(i)+V(i).*Cv(i)-...
%         Cp(i).*(a(1)*W(i)+a(2)*S(i)+Vup(i))+a(7)*Catm(i));