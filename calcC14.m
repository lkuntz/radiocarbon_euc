function [ Cp ] = calcC14( a,x )
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

%Currently neither of these are used, but have been used in older versions
W = x(:,6); %wind anomaly data
T = x(:,7); %Nino3 anomaly Temperatures


Vup_clim = a(1);
Vup_delay = a(2);
Vup_anom_scale = a(3);
Vup_offset = a(4);

Vad_clim = a(5);
Vad_delay = a(6);
Vad_anom_scale = a(7);
Vad_offset = a(8);

Rd_clim = a(9);
Rd_delay = a(10);
Rd_anom_scale = a(11);
Rd_offset = a(12);

H0 = a(13); %initial MLD

ATM_scale = a(14);

months = [0:length(W)-1];
%upwelling velocity
Vup = Vup_clim*sin((months(:)-Vup_delay)*pi/6)-Vup_anom_scale*W+Vup_offset;
%advection velocity
Vad = Vad_clim*sin((months(:)-Vad_delay)*pi/6)-Vad_anom_scale*W+Vad_offset;

%Fraction of upwelled concentration from deep waters
Rdeep = (Rd_clim*sin((months(:)-Rd_delay)*pi/6)-Rd_anom_scale*T+Rd_offset).^2;

h(1) = H0;

for i=1:length(Cn)-1
    
    h(i+1) = h(i)+Vup(i)-Vad(i);

    DCp = ConcentrationToDelta14(Cp(i),1);
    
    %Flux into east pacific box is a mixture of atmospheric penetration
    %(proportional to the difference in atmosphere and previous C14
    %concentration in east Pacific), upwelled deep water, upwelled north
    %and south pacific water, and meridionally advected water.  By mass
    %continuity, the amount of water added equals the amount removed
    %(amount removed is taken to have C14 of the previous time step)
    F(i) = ATM_scale*(Catm(i)-DCp)+Vup(i)*Cd*Rdeep(i) + Vup(i)*(1-Rdeep(i))*(1/3*Cn(i)+2/3*Cs(i))...
        -Vad(i)*(Cp(i))*h(i);
    
    %Concentration of C14 in the next time step is taken by adding the flux
    Cp(i+1) = (h(i)*Cp(i)+F(i))/h(i+1);
end

F=double(F);

Cp = double(Cp);
Cp = 1e4*Cp(:);

end
