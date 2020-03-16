function [GalModeled, GalDICModeled, GalC14, DIC_surface, time] = radiocarbon_model_simulation(dummy, varargin)
    p = inputParser;
    addRequired(p,'dummy');
    addParameter(p,'DIC_deep',2.2,@isnumeric)
    addParameter(p,'DIC_vt',2.05,@isnumeric)
    addParameter(p,'dC14_prebomb',-95,@isnumeric)
    addParameter(p,'ml_depth',25,@isnumeric)
    addParameter(p,'biouptakerate',8*10^-4,@isnumeric)
    addParameter(p,'EUC_data','EUC_nino',@(x) any(validatestring(x, {'EUC_nino', 'EUC_soda', 'EUC_oras'})))
    addParameter(p,'Mixing_Input','data',@(x) any(validatestring(x,{'data', 'Clim', 'Const'})))
    addParameter(p,'WindStress_Input','data',@(x) any(validatestring(x,{'data', 'Clim', 'Const'})))
    addParameter(p,'plot_calibration',false,@islogical)
    
    parse(p,dummy,varargin{:})
    
    DIC_deep = p.Results.DIC_deep;
    DIC_vt = p.Results.DIC_vt;
    dC14_prebomb = p.Results.dC14_prebomb;
    ml_depth = p.Results.ml_depth;
    biouptakerate = p.Results.biouptakerate;
    EUC_data = p.Results.EUC_data;
    Mixing_Input = p.Results.Mixing_Input;
    WindStress_Input = p.Results.WindStress_Input;
    plot_calibration = p.Results.plot_calibration;
    
    %model timestep (1 month)
    dt = 1/12;

    %Set whether projection goes into the future or not
    future_projection = 0;
    
    %DIC estimates from Chai et al 2002 converted to mol/m^3
    DIC_surface = 2.0;
    DIC_subduct = 2.0;

    Rstd = 1/8267; %standard C14

    %Set Radiocarbon concentration of the deep ocean box
    DeepC14 = Delta14ToConcentration(dC14_prebomb, DIC_deep);


    %---------------LOAD TEMPERATURE------------------------------------------
    %data source
    %http://iridl.ldeo.columbia.edu/SOURCES/.GOSTA/
    %load GOSTA climatological temperature data and save necessary info
    [V] = nc_readPH('gosta.clim.nc');
    lat_gosta = V(1).data; %degrees N
    %load lon data (degrees E) and shift so that it starts at 0 and increase to
    %360 (want to avoid having matrix split in the middle of the Pacific)
    lon_gosta = V(3).data; lon_gosta = circshift(lon_gosta,length(lon_gosta)/2);
    lon_gosta(lon_gosta<0) = lon_gosta(lon_gosta<0)+360;
    sst_clim_gosta = V(4).data; %[lon,lat,month]
    %shift the sst values circularly to match the shifted longitude values
    sst_clim_gosta = [sst_clim_gosta(length(lon_gosta)/2+1:end,:,:); sst_clim_gosta(1:length(lon_gosta)/2,:,:)];

    %load GOSTA temperature anomalies
    [V] = nc_readPH('gosta.sst.nc');
    ssta_gosta = V(4).data; %[lon,lat,month]
    %shift the sst anomalies circularly to match shifted longitude values
    ssta_gosta = [ssta_gosta(length(lon_gosta)/2+1:end,:,:); ssta_gosta(1:length(lon_gosta)/2,:,:)];

    %GOSTA starts in Jan 1856 and ends in Dec of 1995
    time_gosta = [1856+1/24:1/12:1995+11/12+1/24];

    %create matrix of nan to store the absolute temperatures
    sst_gosta = nan*ssta_gosta;

    %loop that iterates through each year and adds the monthly climatological
    %values to the anomalies
    for i = 1:length(time_gosta)/12
        sst_gosta(:,:,(i-1)*12+1:i*12) = ssta_gosta(:,:,(i-1)*12+1:i*12)+...
            sst_clim_gosta(:,:,:);
    end


    %-------------------LOAD GALAPAGOES DATA----------------------------------
    %load the galapagos coral data from equatorial western pacific
    %ftp://ftp.ncdc.noaa.gov/pub/data/paleo/coral/east_pacific/urvina14c.txt
    load('GalapagosCoral.csv'); %[year,c14] 90W,0.5S

    %get the SST record from the galapagos coral location.  Average in lat and
    %longitude
    GalT = sst_gosta(54:55,18:19,:); GalT = nanmean(GalT,2); GalT = nanmean(GalT,1);
    GalT = GalT(:);

    %set the time range to interpolate the coral record onto
    if future_projection==1
        time = [1957+1/24:dt:1995]; %when doing future projections
    else
        time = [1957+1/24:dt:GalapagosCoral(1,1)];
    end

    %interpolate the galapagoscoral C14 data onto a monthly time series
    GalC14 = spline(GalapagosCoral(:,1),GalapagosCoral(:,2),time);
    GalC14 = Delta14ToConcentration(GalC14,DIC_surface);

    %plot the data and the interpolated data
    %figure; hold on;
    %plot(GalapagosCoral(:,1),GalapagosCoral(:,2));
    %plot(time,GalC14,'r');


    Nino34 = squeeze(nanmean(nanmean(sst_gosta(39:48,18:19,:),2),1));
    Nino3 = squeeze(nanmean(nanmean(sst_gosta(43:54,18:19,:),2),1));
    Nino3 = double(Nino3);

    % set temperature record to being at start of coral record
    Nino3 = Nino3(find(time_gosta==time(1)):find(time_gosta==time(end)));
    Nino34 = Nino34(find(time_gosta==time(1)):find(time_gosta==time(end)));
    Nino3_anom = Nino3;
    for i=1:12
        Nino3_anom(i:12:end) = Nino3_anom(i:12:end)-nanmean(Nino3_anom(i:12:end));
    end


    %---------------------LOAD WIND STRESS DATA--------------------------------
    %load the wind data and tau data from:
    %http://apps.ecmwf.int/datasets/data/era20c-moda/
    files = dir('ECMWF Data/*nc');
    wind_lat = ncread(strcat('ECMWF Data/',files(1).name),'latitude');
    wind_lon = ncread(strcat('ECMWF Data/',files(1).name),'longitude');

    %region of interest - equatorial E pacific
    Nlat = 5;
    Slat = -5;
    Wlon = 180;
    Elon = 240;

    %region of interest - South Pacific
    N_SPlat = -25;
    S_SPlat = -35;

    %find indices that bound the region of interest
    Nindex = find(wind_lat==Nlat);
    Sindex = find(wind_lat==Slat);
    Windex = find(wind_lon>=Wlon,1,'first');
    Eindex = find(wind_lon>=Elon,1,'first');

    N_SPindex = find(wind_lat==N_SPlat);
    S_SPindex = find(wind_lat==S_SPlat);

    N_NPindex = find(wind_lat==-S_SPlat);
    S_NPindex = find(wind_lat==-N_SPlat);

    Wind_eq = nan(Eindex-Windex+1,Sindex-Nindex+1,12*length(files));
    tau_eq = Wind_eq;
    tau_sp = Wind_eq;
    tau_np = Wind_eq;

    for i=1:length(files)
        WIND = ncread(strcat('ECMWF Data/',files(i).name),'si10');
        WIND = double(WIND);
        TAU = ncread(strcat('ECMWF Data/',files(i).name),'iews');
        TAU = double(TAU);
        Wind_eq(:,:,(i-1)*12+1:i*12) = WIND(Windex:Eindex,Nindex:Sindex,:);
        tau_eq(:,:,(i-1)*12+1:i*12) = TAU(Windex:Eindex,Nindex:Sindex,:);
        tau_sp(:,:,(i-1)*12+1:i*12) = TAU(Windex:Eindex,N_SPindex:S_SPindex,:);
        tau_np(:,:,(i-1)*12+1:i*12) = TAU(Windex:Eindex,N_NPindex:S_NPindex,:);
    end

    Wind_eq = nanmean(Wind_eq,2); Wind_eq = nanmean(Wind_eq,1); Wind_eq = Wind_eq(:);
    tau_eq = nanmean(tau_eq,2); tau_eq = nanmean(tau_eq,1); tau_eq = tau_eq(:);
    tau_sp = nanmean(tau_sp,2); tau_sp = nanmean(tau_sp,1); tau_sp = tau_sp(:);
    tau_np = nanmean(tau_np,2); tau_np = nanmean(tau_np,1); tau_np = tau_np(:);


    time_wind = [1957+1/24:dt:2001-1/24];

    %interpolate to timesteps to match model timesteps
    Wind_eq = spline(time_wind,Wind_eq,time);
    tau_eq = spline(time_wind,tau_eq,time);
    tau_sp = spline(time_wind,tau_sp,time);
    tau_np = spline(time_wind,tau_np,time);


    %using w_ekman = d/dy(Tau/rho*f)
    Omega = 7.2921E-5;
    rho_water = 1025;
    f = 2*Omega*sin(10*pi/180);
    dy = distance(0,180,20,180,referenceEllipsoid('GRS80','m'));


    wE_SP = (tau_sp-tau_eq)/(dy*rho_water*f);
    wE_NP = (tau_np-tau_eq)/(dy*rho_water*f);

    wE = wE_SP+wE_NP;

    wE(wE<0) = 0;

    %figure; plot(wE*3600*24*365/12/20);


    %-------------LOAD SOUTH PACIFIC CORAL DATA-------------------------------%
    %load the data from Rarotonga Coral record in the S Pacific
    %ftp://ftp.ncdc.noaa.gov/pub/data/paleo/coral/east_pacific/rarotonga_14c.txt
    load('RarotongaCoral.csv'); %[year,c14] in S Pacific 159.49W, 21S

    %For multiple measurements at the same time, take the average of the C14
    %data to eliminate duplicates that mess up the interpolation
    RarotongaCoral(398,2) = (RarotongaCoral(398,2)+RarotongaCoral(399,2))/2;

    RarotongaCoral = [RarotongaCoral(1:398,:); RarotongaCoral(400:end,:)];

    %make sure time increases with matrix elements
    RarotongaCoral = flipud(RarotongaCoral);

    %Find the annual minimum value in C14
    %find the start and end index for each year, then find the minimum value.
    %Create a matrix of minimum values and the corresponding time at which that
    %minimum value was measured (we want minimums to keep track of what the
    %subducted signal would have been)
    Year_Raro = floor(RarotongaCoral(:,1));
    start = Year_Raro(1);
    for i = 1:Year_Raro(end)-Year_Raro(1)+1
        startIndex = find(Year_Raro==start,1,'first');
        endIndex = find(Year_Raro==start,1,'last');
        RaroMin(i) = min(RarotongaCoral(startIndex:endIndex,2));
        a = find(min(RarotongaCoral(startIndex:endIndex,2)));
        time_RaroMin(i) = RarotongaCoral(a+startIndex-1,1);
        start = start+1;
    end

    %create a time series onto which to interpolate the minimum rarotoga values
    %and perform the interpolation
    RaroTime = [floor(time_RaroMin(1))+1/24:1/12:ceil(time_RaroMin(end))];
    RarC14 = spline(time_RaroMin,RaroMin,RaroTime);

    % %plot data and interpolated minimums
    % figure; hold on;
    % plot(RarotongaCoral(:,1),RarotongaCoral(:,2));
    % plot(RaroTime,RarC14,'r');

    %Convert to a concentration and plot the results
    RarC14 = Delta14ToConcentration(RarC14,DIC_subduct);

    %Extend the Raratoga time series back in time to 1925, assuming a mean
    %pre-bomb value of -65 per mil
    RarC14 = [linspace(Delta14ToConcentration(-55,DIC_subduct), min(RarC14(1:24)),length([1925+1/24:1/12:RaroTime(1)-1/12])) RarC14];
    RaroTime = [[1925+1/24:1/12:RaroTime(1)-1/12] RaroTime];

    %------------SET UP NORTH PACIFIC DATA------------------------------------%
    % %load data from the North pacific coral
    % %Druffel (1987) Bomb radiocarbon in the Pacific: Annual and seasonal
    % %timescale variations. Journal of Marine Research, 45, 667-698
    %load('FrenchFrigateCoral.csv'); %[year,c14] in N Pacific 166W,23.7N

    %use the rarotoga coral record to also represent the north pacific assuming
    %the shape of the two curves is the same, just offset in time by 2 years
    %and with a larger magnitude in the northern hemisphere
    FFC14 = (RarC14-RarC14(1))*1.18+RarC14(1); %scale the data, but still have it start at -65 per mil
    FFTime = RaroTime-2;


    %----------CALCULATE C14 Concentration in the EUC-------------------------%
    %find start time and and end time that align with Galapagoes Record
    RarIndex = find(RaroTime==time(1));
    FFIndex = find(FFTime == time(1));
    RFIndex_start = find(FFTime == RaroTime(1));
    RFIndex_end = find(RaroTime==FFTime(end));

    %Calculate EUC C14 as smooth signal of N and S Pacific (20% and 80%
    %respectively)
    EUCcon = .2*FFC14(RFIndex_start:end)+.8*RarC14(1:RFIndex_end);
    EUCcon = tsmovavg(EUCcon,'s',12*12);

    %Delay due to transport of signal
    transit =24;
    EUCc14 = EUCcon(FFIndex-transit:FFIndex-transit+length(GalC14)-1);

    %-----------------LOAD ATMOSPHERIC C14 DATA-------------------------------%
    %load data of atmospheric record of C14 in the Northern Hemisphere
    %http://cdiac.ess-dive.lbl.gov/ftp/trends/co2/vermunt.c14
    %atm radiocarbon measure at Vermunt, Austria
    load('VermuntC14.csv'); %[month,year,d13c,Dc14]
    time_verATM = 1900+VermuntC14(:,2)+(VermuntC14(:,1)-.5)/12;
    verC14 = VermuntC14(:,4);

    %take average when multiple measurements for one month
    [C,ia,idx] = unique(time_verATM,'stable');
    verC14 = accumarray(idx,verC14,[],@mean); 
    verc13 = accumarray(idx,VermuntC14(:,3),[],@mean);
    time_verC14 = C;

    %load data of atmospheric record of C14 in the Southern Hemisphere
    %http://cdiac.esd.ornl.gov/ftp/trends/co2/welling.195
    %atm radiocarbon measure at Wellington, New Zealand
    load('welling.195.csv'); %[year,month,c13,c14]

    %time at which measurements were made
    timeAtm = welling_195(:,1)+1900+(welling_195(:,2)-.5)/12;

    %take average when multiple measurements for one month
    [C,ia,idx] = unique(timeAtm,'stable');
    atmC14 = accumarray(idx,welling_195(:,4),[],@mean); 
    atmc13 = accumarray(idx,welling_195(:,3),[],@mean);
    timeAtm = C;


    %average in the data from the N hemisphere
    for i=1:length(verC14)
        ind = find(timeAtm==time_verC14(i));
        atmC14(ind) = atmC14(ind)/2+verC14(i)/2;
        atmc13(ind) = atmc13(ind)/2+verc13(i)/2;
    end

    %interpolate the atmospheric measurements onto a monthly time series (note
    %that the annual cycle of C14 changes has been removed from this)
    AtmC14 = spline(timeAtm,atmC14,time);
    AtmC13 = spline(timeAtm,atmc13,time);

    %plot the data
    % figure; plot(time,AtmC14);

    %--------------------LOAD ATMOSPHERIC CO2 DATA----------------------------%
    %From: ftp://aftp.cmdl.noaa.gov/products/trends/co2/co2_mm_mlo.txt
    load('CO2_ppm_Hawaii.csv'); %[year,month,time,ppm]
    time_CO2 = CO2_ppm_Hawaii(:,3);
    CO2_Hawaii = CO2_ppm_Hawaii(:,4);
    time_CO2 = [1957;time_CO2];
    CO2_Hawaii = [mean(CO2_Hawaii(1:12));CO2_Hawaii];
    AtmCO2 = interp1(time_CO2,CO2_Hawaii,time);
    % figure; plot(time,AtmCO2);

    %--------------------CALCULATE MAX VELOCITIES-----------------------------%
    %Determine which dataset to use to get maximum velocity of EUC

    if strcmp(EUC_data,'EUC_nino')
        vEUC = (Nino34'-31.01)/-.03678; %unsmoothed, straight linear fit from 
        %CESM - Temperature based Nino34 Regression
    else
        if strcmp(EUC_data,'EUC_soda')
            %data from SODA reanalysis
            %http://iridl.ldeo.columbia.edu/SOURCES/.CARTON-GIESE/.SODA/.v2p2p4/.u/
            uvel = ncread('SODA_1957-2008.nc','u');
            ulat = ncread('SODA_1957-2008.nc','lat');
            ulon = ncread('SODA_1957-2008.nc','lon');
        else strcmp(EUC_data,'EUC_oras')
            %data from ORAS reanalysis
            uvel = ncread('/Users/Lauren/Documents/Harvard EPS/Sulfate Aerosols/Pacific Simulation/ORAS4-ZonalV.nc','UO');
            time_uvel = [1957+8/12:1/12:2015+7/12];
            ulat = ncread('/Users/Lauren/Documents/Harvard EPS/Sulfate Aerosols/Pacific Simulation/ORAS4-ZonalV.nc','LAT89_92');
            ulon = ncread('/Users/Lauren/Documents/Harvard EPS/Sulfate Aerosols/Pacific Simulation/ORAS4-ZonalV.nc','LON120_291');
            uvel = cat(4,nanmean(uvel(:,:,:,5:12:end),4),nanmean(uvel(:,:,:,6:12:end),4),...
                nanmean(uvel(:,:,:,7:12:end),4),nanmean(uvel(:,:,:,8:12:end),4),...
                nanmean(uvel(:,:,:,9:12:end),4),nanmean(uvel(:,:,:,10:12:end),4),...
                nanmean(uvel(:,:,:,11:12:end),4),nanmean(uvel(:,:,:,12:12:end),4),...
                uvel); %extend record back using climatology
        end
        %create grid of lat longs for interpolation
        [uulat, uulon] = meshgrid(ulat,ulon);
        uulat = double(uulat);
        uulon = double(uulon);
        uvel = double(uvel);

        %interpolate data to 220E, 0N
        for i=1:size(uvel,4)
            for j=1:size(uvel,3)
                tointerp = squeeze(uvel(:,:,j,i));
                interpolant = scatteredInterpolant(uulon(:),uulat(:),tointerp(:));
                u_220(j,i) = interpolant(220,0);
            end
        end
        %find max value at each timestep
        for i=1:size(u_220,2)
            vEUC(i) = max(u_220(:,i));
        end
        vEUC = vEUC*100; %convert to cm/s
    end

    %calculate ratio of mixing that occurs in thermocline box
    rDeep = vEUC.^2*.2/mean(vEUC(1:length(GalC14)).^2);

    %calculate Ekman Pumping per month as fraction of ml_depth
    vReplace = 1*wE*3600*24*365/12/ml_depth;
    vReplace(vReplace<0) = 0;

    %Initial values for the model
    Gal0 = GalC14(1);
    GalModeled = Gal0;
    GalDICModeled = DIC_surface;

    %Calculate climatologies of upwelling and mixing
    vReplaceClim = [];
    rDeepClim = [];
    for i=1:12
        vReplaceClim(i) = nanmean(vReplace(i:12:end));
        rDeepClim(i) = nanmean(rDeep(i:12:end));
    end

    %if using constant or climatological mixing, change rDeep
    if strcmp(Mixing_Input,'Const')
        rDeep = nanmean(rDeepClim)*ones(1,length(time));
    elseif strcmp(Mixing_Input,'Clim')
        rDeep = repmat(rDeepClim,[1 50]);
    end

    %if using constant or climatological upwelling change vReplace
    if strcmp(WindStress_Input,'Const')
        vReplace = nanmean(vReplaceClim)*ones(1,length(time));
    elseif strcmp(WindStress_Input,'Clim')
        vReplace = repmat(vReplaceClim,[1 50]);
    end


    %----------------------MODEL SIMULATION-----------------------------------%

    %monthly timestep through the model
    for i=1:length(Nino3)
        %parameter calculation for DIC solubility - air-gas exchange
        [Ko,k] = DIC_solubility(Nino3(i),34.78,Wind_eq(i));

        %invasion of CO2 and C14 from atmosphere
        InvasionCO2 = k*Ko*AtmCO2(i)/1e6; %[mol/m^2/month]
        InvasionC14 = InvasionCO2*(Delta14Todelta(AtmC14(i),AtmC13(i))/1000+1)*Rstd;

        %biological update of DIC and C14
        BioUptake = biouptakerate*365/12; %From Chai et al Fig. 7
        BioUptakeC14 = BioUptake*GalModeled(i)/GalDICModeled(i);

        %outgasing of CO2 and C14 from ocean
        EvasionCO2 = k*Ko*carbonate_calculation(Nino3(i),GalDICModeled(i))/1e6;
        EvasionC14 = EvasionCO2*GalModeled(i)/GalDICModeled(i);

        %concentration of DIC and C14 in VT box
        upC14(i) = (rDeep(i)*DeepC14+rDeep(i)*GalModeled(i)+EUCc14(i)*(1-2*rDeep(i)));
        upDIC(i) = (rDeep(i)*DIC_deep+rDeep(i)*GalDICModeled(i)+DIC_vt*(1-2*rDeep(i)));

        GalModeled(i+1) = GalModeled(i)*(1-vReplace(i)-rDeep(i))+upC14(i)*(vReplace(i)+rDeep(i))+(InvasionC14-EvasionC14)/ml_depth-BioUptakeC14;
        GalDICModeled(i+1) = GalDICModeled(i)*(1-vReplace(i)-rDeep(i))+upDIC(i)*(vReplace(i)+rDeep(i))+(InvasionCO2-EvasionCO2)/ml_depth-BioUptake;
    end
    
    if plot_calibration

        %---------------MODEL CALIBRATION TO PRE-BOMB VALUES----------------------%
        vReplace = repmat(vReplaceClim,[1 15]); %use climatology of upwelling
        rDeep = repmat(rDeepClim,[1 15]); %use climatology of mixing with deep

        %Set initial conditions for model
        Gal0 = Delta14ToConcentration(-72, DIC_surface);
        GalModeled = Gal0;
        GalDICModeled = DIC_surface;

        %monthly timestep through model
        for i=1:15*12
            [Ko,k] = DIC_solubility(Nino3(1),34.78,Wind_eq(1));

            InvasionCO2 = k*Ko*AtmCO2(1)/1e6; %[mol/m^2/month]
            InvasionC14 = InvasionCO2*(Delta14Todelta(AtmC14(1),AtmC13(1))/1000+1)*Rstd;

            BioUptake = 20.09*10^-3*365/12/ml_depth; %From Chai et al Fig. 7
            BioUptakeC14 = BioUptake*GalModeled(i)/GalDICModeled(i);


            EvasionCO2 = k*Ko*carbonate_calculation(Nino3(i),GalDICModeled(i))/1e6;
            EvasionC14 = EvasionCO2*GalModeled(i)/GalDICModeled(i);

            upC14(i) = (rDeep(i)*DeepC14+rDeep(i)*GalModeled(i)+EUCc14(1)*(1-2*rDeep(i)));
            upDIC(i) = (rDeep(i)*DIC_deep+rDeep(i)*GalDICModeled(i)+DIC_vt*(1-2*rDeep(i)));

            GalModeled(i+1) = GalModeled(i)*(1-vReplace(i)-rDeep(i))+upC14(i)*(vReplace(i)+rDeep(i))+(InvasionC14-EvasionC14)/ml_depth-BioUptakeC14;
            GalDICModeled(i+1) = GalDICModeled(i)*(1-vReplace(i)-rDeep(i))+upDIC(i)*(vReplace(i)+rDeep(i))+(InvasionCO2-EvasionCO2)/ml_depth-BioUptake;
        end

        f = figure;
        set(f,'Units','normalized');
        set(f,'Position',[0 0 1 1]);
        plot([-1:1/12:14-1/12],ConcentrationToDelta14(GalModeled(1:end-1),GalDICModeled(1:end-1)),'LineWidth',2);
        xlabel('Year'); ylabel('\Delta^{14}C');
        set(gca,'FontSize',16);
        xlim([0 10]);

        f = figure;
        set(f,'Units','normalized');
        set(f,'Position',[0 0 1 1]);
        plot([-1:1/12:14-1/12],GalDICModeled(1:end-1),'LineWidth',2);
        xlabel('Year'); ylabel('DIC Concentration [mol/m^3]');
        set(gca,'FontSize',16);
        xlim([0 10]);
    end
end
