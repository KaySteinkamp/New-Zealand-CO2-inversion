function data = read_data(invcfg, obscfg, sources)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Reads observational data from BHD, LAU and RBM stations, and
% puts together a datafile for a given period usable in the inversion.
% 
% It returns the names of stations, observations, variances etc. 
% Returned quantities are stored in the structure "data".
% 
% Author: Kay Steinkamp
% Date: Jan 2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set adjustable parameter
if isfield(obscfg,'scale')
    adjust = obscfg.scale;     % optimized value
else
    adjust = obscfg.iniscale;  % initial (default) value
end


% open and read the station data files
% ----------------------------------------------------------------------
if obscfg.syntheticObs
    
    % Load synthetic CO2 anomalies (daily @15:00 LT and @13:00 LT).
    % NOTE, anomalies include response to land/ocean/FF flux
    load([invcfg.workdir obscfg.Fsynobs],'synobs13h','synobs15h')
        
    switch invcfg.type
        case 'monthly'
            Xco2 = 'xfrac';
        case 'weekly'
            Xco2 = 'xfrac_wi';
        otherwise
            error('Use "monthly" or "weekly" as inversion type')
    end
    
    % Weighted N/S baseline differs among sites, but the old STL does not
    if obscfg.useWtBsl
        Xbsl = 'xbslwt';
        nreplicate = 1;
    else
        Xbsl = 'xstl';
        nreplicate = 3;
    end
    %%% TEST %%%
%     Xbsl = 'xstl';
%     nreplicate = 3;
    %%%%%%%%%%%%
    pX = [];
    stl = [];
    
    if obscfg.include13h
        % Daily synthetic CO2 @13:00 LT, including baseline and
        % fossil response, in ppm
        if ~isequal(synobs13h.site,{'BHD'  'LAU'  'RBM'})
            error('Unexpected order of stations')
        end
        if obscfg.synthDiCy
            error('No synthetic data available for DiCY case at 13h')
        end
        if isempty(pX)
            stl = repmat(synobs13h.(Xbsl),[1 nreplicate]); % baseline
            pX = synobs13h.(Xco2) + stl; % BHD, LAU, RBM
        else
            stl = [stl; repmat(synobs13h.(Xbsl),[1 nreplicate])];
            pX = [pX; synobs13h.(Xco2) + repmat(synobs13h.(Xbsl),[1 nreplicate])];
        end
    end
    
    if obscfg.include15h
        % Daily synthetic CO2 @15:00 LT, including baseline and
        % fossil response, in ppm
        if ~isequal(synobs15h.site,{'BHD'  'LAU'  'RBM'})
            error('Unexpected order of stations')
        end
        if obscfg.synthDiCy
            % This is the test case with diurnal cycle
            stl = repmat(synobs15h.(Xbsl),[1 nreplicate]);
            %pX = synobs15h.xfracCtrl_wi + stl; % BHD, LAU, RBM
            %pX = synobs15h.xfracCtrlLa + stl; % BHD, LAU, RBM
            pX = synobs15h.xfracLa + stl; % BHD, LAU, RBM
        else
            if isempty(pX)
                stl = repmat(synobs15h.(Xbsl),[1 nreplicate]); % baseline
                pX = synobs15h.(Xco2) + stl; % BHD, LAU, RBM
            else
                stl = [stl; repmat(synobs15h.(Xbsl),[1 nreplicate])];
                pX = [pX; synobs15h.(Xco2) + repmat(synobs15h.(Xbsl),[1 nreplicate])];
            end
        end
    end
    
    % Assign uncertainty
%     pXunc = repmat(obscfg.synobsUnc,size(pX));
%     stlunc = repmat(obscfg.synobsUnc,size(stl));
    if obscfg.synthDiCy
        pXunc = synobs15h.xfracUnc;
        stlunc = synobs15h.xbslwtUnc;
    else
        pXunc = repmat(obscfg.synobsUnc,size(pX));
        stlunc = repmat(obscfg.synobsUnc,size(stl));
    end
    
    % Metadata
    sitenm = synobs15h.site; 
    day2mon = synobs15h.day2mon;
    day = datenum(synobs15h.day);
    isDS = zeros(size(day));
    
else
    
    % Load data set and - if needed - weighted baseline
    load([invcfg.workdir obscfg.Fobs],'Data3sites','Header3sites')
    if obscfg.useWtBsl
        load([invcfg.workdir obscfg.FwtBsl],'bsl13h','bsl15h')
    end
    
    isDS = [];
    day = [];
    stl = [];
    pX = [];
    pXunc = [];
    
    if obscfg.include13h
        % Daily CO2 @13:00 LT and uncertainty, in ppm
        clear *_tmp
        pick = find((Data3sites(:,strcmp(Header3sites,'StartHour(NZST)'))==13 &...
                     Data3sites(:,strcmp(Header3sites,'Daylight Saving'))==0) ...
                  | (Data3sites(:,strcmp(Header3sites,'StartHour(NZST)'))==14 &...
                     Data3sites(:,strcmp(Header3sites,'Daylight Saving'))==1)); %#ok<NODEF>
        isDS_tmp = Data3sites(pick,strcmp(Header3sites,'Daylight Saving'));
        day_tmp = Data3sites(pick,strcmp(Header3sites,'numTime'));
        if obscfg.useWtBsl
            % weighted N/S baseline
            stl_tmp = bsl13h.XbslWeight;
            stlunc_tmp = bsl13h.XbslWeightSD;
        else
            % BHD baseline
            stl_tmp = Data3sites(pick,strcmp(Header3sites,'BHD: STLhourly')); 
            stlunc_tmp = Data3sites(pick,strcmp(Header3sites,'BHD: STLhourly-StDev')); 
        end
        pX_tmp(:,1) = Data3sites(pick,strcmp(Header3sites,'BHD: co2')); % BHD
        pXunc_tmp(:,1) = Data3sites(pick,strcmp(Header3sites,'BHD: sd'));
        pX_tmp(:,2) = Data3sites(pick,strcmp(Header3sites,'LAU: CO2_1')); % LAU
        pXunc_tmp(:,2) = Data3sites(pick,strcmp(Header3sites,'LAU: Sigma_1'));
        pX_tmp(:,3) = Data3sites(pick,strcmp(Header3sites,'RBM: CO2_Mean')); % RBM
        % pX_tmp(:,3) = Data3sites(pick,strcmp(Header3sites,'RBM: CO2_1')); 
        pXunc_tmp(:,3) = Data3sites(pick,strcmp(Header3sites,'RBM: Sigma_1'));
        
        % exclude flagged RBM data
        flagged = Data3sites(pick,strcmp(Header3sites,'RBM: flag'));
        pX_tmp(flagged==1,3) = nan;
        pXunc_tmp(flagged==1,3) = nan;
        
        if isempty(pX)
            isDS = isDS_tmp;
            day = day_tmp;
            stl = stl_tmp;
            stlunc = stlunc_tmp;
            pX = pX_tmp;
            pXunc = pXunc_tmp;
        else
            stl = [stl; stl_tmp];
            stlunc = [stlunc; stlunc_tmp];
            pX = [pX; pX_tmp];
            pXunc = [pXunc; pXunc_tmp];
        end
    end
    
    if obscfg.include15h
        % Daily CO2 @15:00 LT and uncertainty, in ppm
        clear *_tmp
        pick = find((Data3sites(:,strcmp(Header3sites,'StartHour(NZST)'))==15 &...
                     Data3sites(:,strcmp(Header3sites,'Daylight Saving'))==0) ...
                  | (Data3sites(:,strcmp(Header3sites,'StartHour(NZST)'))==16 &...
                     Data3sites(:,strcmp(Header3sites,'Daylight Saving'))==1)); 
        isDS_tmp = Data3sites(pick,strcmp(Header3sites,'Daylight Saving'));
        day_tmp = Data3sites(pick,strcmp(Header3sites,'numTime'));
        if obscfg.useWtBsl
            % weighted N/S baseline
            stl_tmp = bsl15h.XbslWeight;
            stlunc_tmp = bsl15h.XbslWeightSD;
        else
            % BHD baseline
            stl_tmp = Data3sites(pick,strcmp(Header3sites,'BHD: STLhourly'));
            stlunc_tmp = Data3sites(pick,strcmp(Header3sites,'BHD: STLhourly-StDev'));
        end
        pX_tmp(:,1) = Data3sites(pick,strcmp(Header3sites,'BHD: co2')); % BHD
        pXunc_tmp(:,1) = Data3sites(pick,strcmp(Header3sites,'BHD: sd'));
        pX_tmp(:,2) = Data3sites(pick,strcmp(Header3sites,'LAU: CO2_1')); % LAU
        pXunc_tmp(:,2) = Data3sites(pick,strcmp(Header3sites,'LAU: Sigma_1'));
        pX_tmp(:,3) = Data3sites(pick,strcmp(Header3sites,'RBM: CO2_Mean')); % RBM
        % pX_tmp(:,3) = Data3sites(pick,strcmp(Header3sites,'RBM: CO2_1')); 
        pXunc_tmp(:,3) = Data3sites(pick,strcmp(Header3sites,'RBM: Sigma_1'));
        
        % exclude flagged RBM data
        flagged = Data3sites(pick,strcmp(Header3sites,'RBM: flag'));
        pX_tmp(flagged==1,3) = nan;
        pXunc_tmp(flagged==1,3) = nan;
        
        if isempty(pX)
            isDS = isDS_tmp;
            day = day_tmp;
            stl = stl_tmp;
            stlunc = stlunc_tmp;
            pX = pX_tmp;
            pXunc = pXunc_tmp;
        else
            stl = [stl; stl_tmp];
            stlunc = [stlunc; stlunc_tmp];
            pX = [pX; pX_tmp];
            pXunc = [pXunc; pXunc_tmp];
        end
    end
    
    % Metadata
    sitenm = {'BHD','LAU','RBM'};
    tmp = datevec(day);
    day = datenum(tmp(:,1:3)); % keep only day - no times
    day2mon = [];
    for y = invcfg.period(1):invcfg.period(2)
        day2mon = [day2mon; tmp(tmp(:,1)==y,2)+(y-invcfg.period(1))*12];
    end

end


% Adjust data & uncertainty
% ----------------------------------------------------------------------
%%% TESTs
% pXunc = 2.5*ones(size(pX));
% pXunc = pXunc .* 6;
% % randomly omit data:
% exclude = randperm(numel(pX),round(0.7*numel(pX)));
% pX(exclude) = nan;
% pXunc(exclude) = nan;
% pX(:,2) = pX(:,2) + (0.5+1.4*(1:-1/730:0))'; % counter too high baseline
%%%

% Adjust CO2 to NOAA CCL, based on Oct 2014 Aniwaniwa tank suite
if obscfg.adjustCCL
    bhd2noaa = @(x) 1.00286417.*x - 1.15254309;
    lau2noaa = @(x) 0.99904604.*x + 0.45969191;
    rbm2noaa = @(x) 0.997412.*x + 1.160247 ;
else
    bhd2noaa = @(x) x;
    lau2noaa = @(x) x;
    rbm2noaa = @(x) x;
end
pX(:,1) = bhd2noaa(pX(:,1));
pX(:,2) = lau2noaa(pX(:,2));
pX(:,3) = rbm2noaa(pX(:,3));
 
% exclude/include stations
if ~obscfg.includeBHD
    pX(:,1) = nan;
    pXunc(:,1) = nan;
end
if ~obscfg.includeLAU
    pX(:,2) = nan;
    pXunc(:,2) = nan;
end
if ~obscfg.includeRBM
    pX(:,3) = nan;
    pXunc(:,3) = nan;
end


% Apply uniform noise and bias if configured
% ----------------------------------------------------------------------
pX = pX + obscfg.biasX + obscfg.uniNoiseX*(rand(size(pX))/0.5 - 1);


% Subtract baseline and combine obs & baseline uncertainty
% Also apply baseline bias if configured
% ----------------------------------------------------------------------
if obscfg.useWtBsl
    stlBias = obscfg.BslBiasSig*stlunc;
    pX = pX - (stl + stlBias);
    pXunc = sqrt(pXunc.^2 + stlunc.^2);
else
    pX = pX - stl;
end


% Apply minimum uncertainty and adjust for chi^2=1
% ----------------------------------------------------------------------
pXunc(pXunc < obscfg.minUnc) = obscfg.minUnc;
%pXunc = sqrt(pXunc.^2 + obscfg.minUnc.^2);

% 'adjust' will differ from one if chi^2 optimisation is on
pXunc = pXunc .* adjust;


% Exclude data outside 3 standard deviations
% ----------------------------------------------------------------------
% % stdlim = 3*nanstd(pX);
% for s = 1:3
%     stdlim = 3*std(pX(~isnan(pX(:,s)),s));
%     pX(abs(pX(:,s))>stdlim,s) = nan;
% end


% Subtract fossil fuel response
% ----------------------------------------------------------------------    
load([invcfg.workdir obscfg.FrespFoss],'fresp13h','fresp15h')

% Ensure consistent day and site arrays
if ~obscfg.syntheticObs && obscfg.include13h && ...
        any(abs(day - (datenum(fresp13h.day))) > 1e-9)
    error('Inconsistent day arrays')
end
if ~obscfg.syntheticObs && obscfg.include15h && ...
        any(abs(day - (datenum(fresp15h.day))) > 1e-9)
    error('Inconsistent day arrays')
end
if ~isequal(fresp13h.site,{'BHD','LAU','RBM'}) ||...
        ~isequal(fresp15h.site,{'BHD','LAU','RBM'})
    error('Unexpected order of stations')
end

xFF = [];
if obscfg.presubFoss
    
    if obscfg.include13h
        xFF = [xFF; fresp13h.xfracL];
    end
    
    if obscfg.include15h
        xFF = [xFF; fresp15h.xfrac];
    end
    
    pX = pX - xFF; 
end


% Smooth timeseries using 7 day span
% (to match up timescale of data/sources=weekly/weekly)
% ----------------------------------------------------------------------
if obscfg.smoothObs
    nobsperday = size(pX,1)/length(day);
    if nobsperday == 1
        for i = 1:3
            pX(:,i) = smooth(pX(:,i),3,'moving');
            pXunc(:,i) = smooth(pXunc(:,i),3,'moving');
        end
    elseif nobsperday == 2
        nd = length(day);
        for i = 1:3
            pX(1:nd,i) = smooth(pX(1:nd,i),3,'moving');
            pX(nd+1:end,i) = smooth(pX(nd+1:end,i),3,'moving');
            pXunc(1:nd,i) = smooth(pXunc(1:nd,i),3,'moving');
            pXunc(nd+1:end,i) = smooth(pXunc(nd+1:end,i),3,'moving');
        end
    else
        error('Unexpexted number of observations per day')
    end
end


% Convert ppm to g CO2 m-3, using NAME-NZLAM temperature and pressure
% ----------------------------------------------------------------------
T = [];
P = [];

% @13-14h
if obscfg.include13h
    if strcmp(fresp13h.tempUnit,'deg C') && strcmp(fresp13h.presUnit,'Pa')
        T = [T; fresp13h.tempL + 273.15]; % K
        P = [P; fresp13h.presL]; % Pa
    else
        error('Unexpected units for model temperature and/or pressure')
    end
end

% @15-16h
if obscfg.include15h
    if strcmp(fresp15h.tempUnit,'deg C') && strcmp(fresp15h.presUnit,'Pa')
        T = [T; fresp15h.temp + 273.15]; % K
        P = [P; fresp15h.pres]; % Pa
    else
        error('Unexpected units for model temperature and/or pressure')
    end
end


pC = (44/8.31) * (pX*1e-6) .* P ./ T;
pCunc = (44/8.31) * (pXunc*1e-6) .* P ./ T;


% Assign large uncertainty for days/sites with data, but without sd
% ----------------------------------------------------------------------
nosd = ( isnan(pXunc) & ~isnan(pX) );
if ~isequal(nosd,(isnan(pCunc)&~isnan(pC)))
    error('Something is really wrong here')
end
pXunc(nosd) = 10; % 10 ppm ~ 0.02 g CO2 m-3
pCunc(nosd) = (44/8.31) * (pXunc(nosd)*1e-6) .* P(nosd) ./ T(nosd);


% Reshape day*site matrix to vector form
% ----------------------------------------------------------------------
nsites = length(sitenm);
nday = length(day);
ndat = numel(pC);
ntimes = ndat/nsites/nday;
dayV = repmat(day,[1 nsites ntimes]);
dayV = reshape(dayV,ndat,1);
pCV = reshape(pC,ndat,1);
pCuncV = reshape(pCunc,ndat,1);
pXV = reshape(pX,ndat,1);
pXuncV = reshape(pXunc,ndat,1);
TV = reshape(T,ndat,1);
PV = reshape(P,ndat,1);


% Exclude days/sites with unavailable data (only vector forms)
% ----------------------------------------------------------------------
iavail = find(~isnan(pCV));
if ~isequal(iavail,find(~isnan(pXV)))
    error('Something is really wrong here')
end
pCV = pCV(iavail);
pXV = pXV(iavail);
pCuncV = pCuncV(iavail);
pXuncV = pXuncV(iavail);

TV = TV(iavail);
PV = PV(iavail);


% Temporal smoothing of fluxes as extra 'observations'
% ----------------------------------------------------------------------
w2wdev_land = invcfg.smoothWland; % kg CO2 m-2 yr-1
w2wdev_oce = invcfg.smoothWoce; % kg CO2 m-2 yr-1

if (w2wdev_land>0 || w2wdev_oce>0) && strcmp(invcfg.type,'weekly')
    
    % number of regions x weeks
    nweek = 48*(diff(invcfg.period)+1);
    nextraObs = 25*nweek;
    
    % there are nweek (=144) weeks for each region
    load([invcfg.workdir 'NAMEgrid.mat'],'NAMEgrid')
    w2wdev_reg(1:19) = w2wdev_land*NAMEgrid.regarea(1:19)*1e-9; % Tg CO2 yr-1
    w2wdev_reg(20:25) = w2wdev_oce*NAMEgrid.regarea(20:25)*1e-9; % Tg CO2 yr-1
    w2wdev = nan(nextraObs,1);
    for r = 1:25
        w2wdev((r-1)*nweek+1:r*nweek) = w2wdev_reg(r);
    end
    if any(isnan(w2wdev))
        error('failed to set up temporal smoothing')
    end

    % apply Gaussian smoother with fixed week-to-week stdev (Tg CO2 yr-1)
    % add equation for each week and each region
    pXV = [pXV; zeros(nextraObs,1)];
    pXuncV = [pXuncV; w2wdev];
    pCV = [pCV; zeros(nextraObs,1)];
    pCuncV = [pCuncV; w2wdev];
    
    TV = [TV; nan(nextraObs,1)];
    PV = [PV; nan(nextraObs,1)];
    
else
    nextraObs = 0;
end


% Annual constraint for Open Ocean (OO) and Australia (Oz)
% ----------------------------------------------------------------------
annunc_pc = invcfg.OzOOAnn_pc; % per cent

if (annunc_pc > 0) && invcfg.oceanprior
    
    % number of years x Oz and OO regions ([16 20:25])
    Iann = [16 20:25];
    nann = length(Iann)*(diff(invcfg.period)+1);
    nextraObs = nextraObs + nann;
    
    % number of sources per region and year
    switch invcfg.type
        case 'weekly'
            ns = 48;
        case 'monthly'
            ns = 12;
    end
    
    % write annual uncs at correct locations in array
    annsrc = []; annunc = [];
    for r = Iann
        for y = 1:diff(invcfg.period)+1
            asrc = mean(sources.prior((y-1)*ns+1:y*ns,r)); % Tg CO2 yr-1
            if r == 16
                annunc = [annunc; 1000]; % Tg CO2 yr-1
            else
                annunc = [annunc; abs(asrc)*annunc_pc/100]; % Tg CO2 yr-1
            end
            annsrc = [annsrc; asrc];
        end
    end

    % add equation for each year and each region
    pXV = [pXV; ns*annsrc];
    pXuncV = [pXuncV; ns*annunc];
    pCV = [pCV; ns*annsrc];
    pCuncV = [pCuncV; ns*annunc];
    
    TV = [TV; nan(nann,1)];
    PV = [PV; nan(nann,1)];

end


% Store in data structure
% ----------------------------------------------------------------------
data.nsites = nsites;
data.sitename = sitenm; 
data.sitecoord = [-41.4083 -45.03 -38.3193;
                  174.8710 169.68 176.3882];     
data.ntimes = ntimes;
data.ndat = ndat;
data.day = day;
data.dayV = dayV;
data.daylightSaving = isDS;
data.Note = 'Dates given in NZST (daylight saving is accounted for)';
data.day2mon = day2mon;
data.baselineX = stl;
data.baselineXunc = stlunc;
data.FFrespX = xFF;
data.nextra = nextraObs;
data.neqns = length(iavail) + nextraObs;
data.iavail = iavail;
data.type = 'Daily anomalies';
if obscfg.include13h
    data.type = [data.type ' @13:00 LT'];
end
if obscfg.include15h
    data.type = [data.type ' @15:00 LT'];
end
data.priorX = pX; 
data.priorXunc = pXunc;
data.priorXV = pXV; 
data.priorXuncV = pXuncV;
data.unitsX = 'ppm';
data.T = T;
data.TV = TV;
data.unitsT = 'K';
data.P = P;
data.PV = PV;
data.unitsP = 'Pa';
data.priorC = pC; 
data.priorCunc = pCunc;
data.priorCV = pCV; 
data.priorCuncV = pCuncV;
data.unitsC = 'g CO2 m-3';
data.chi2adjust = adjust;


% Display some information
% ----------------------------------------------------------------------
if invcfg.chatty
    fprintf('\nStation datasets have been read for %g-%g\n',...
        invcfg.period(1),invcfg.period(2));
    fprintf('Number of stations             : %g\n',data.nsites);
    fprintf('Number of days in period       : %g\n',nday);
    fprintf('Number of data per day         : %g\n',ntimes);
    fprintf('Number of available data points: %g\n',data.neqns);
    fprintf('Range of data (anomalies)      : %5.2f to %5.2f ppm\n',...
        min(pX(:)),max(pX(:)));
    fprintf('Range of data uncertainty      : %5.2f to %5.2f ppm\n',...
        min(pXunc(:)),max(pXunc(:)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%