function frech = makefrech(invcfg, obscfg, nsites, ntimes, ndat, ...
                           iavail, sources)   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Builds the equation part of the Frechet/Jacobi matrix for the NZ
% inversion. Additional a priori constraints are added as well.
%
% First, the Green's functions (i.e. the model responses at the
% observation sites) are built.
% From that the Frechet (i.e. the Jacobian) is built: 
% each matrix element is just the Green's function for the particular 
% data location multiplied by a "recipe" for converting a particular
% source type into a particular data type. In the case of the NZ
% inversion, those scaling factors simply equal one (only daily data  
% and monthly sources are used).
%
% Author: Kay Steinkamp
% Date: Jan 2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global weekspan Nday dayrange

nreg = sources.nreg;
nvar = sources.nvar;
nsrc = nvar - nsites; 
nmon = sources.nmon;
nweek = sources.nweek;


% Read daily responses at sites/locations (@13:00 LT and/or @15:00 LT)
% ----------------------------------------------------------------------
% Load Green's Functions (in g CO2 m-3) for each day in inv period (more 
% (more precisely, each 13:00-14:00 and/or 15:00-16:00 LT measurement), 
% and for each site.
% ----------------------------------------------------------------------
load([invcfg.workdir invcfg.FGreensFun],'resp13h','resp15h');

% get release days 
% (should be 5 for each obs day: obs day plus 4 days of backtrajectories)
relDays = resp15h.relDays;
if ~isequal(relDays,resp13h.relDays)
    error('Inconsistent release history')
end 


% Make consistency checks and define day-to-month and day-to-week vector
% ----------------------------------------------------------------------
if ~isequal(resp13h.day2mon,resp15h.day2mon) || ...
        ~isequal(resp13h.day,resp15h.day)
    error('Inconsistent month/day arrays')
end

nday = length(resp13h.day);
if ~isequal(Nday,nday)
    error('Inconsistent number of days')
end
firstDay = datenum(dayrange{1});

if ndat ~= nday*nsites*ntimes
    error('Confused by dimensions!')
end

% define vector translating day to month
day2mon = resp15h.day2mon;
if ~isequal(nmon,max(day2mon))
    error('Inconsistent number of months')
end 

% define vector translating day to week
if strcmp(invcfg.type,'weekly')
    if ~isequal(nweek,4*nmon)
        error('Inconsistent number of weeks')
    end
    day2week = zeros(nweek,1);
    for d = 1:nday
        day2week(d) = find(d>weekspan(1,:) & d<=weekspan(2,:));
    end
else
    day2week = [];
end 
        
        
% Choose within-region source pattern
% ----------------------------------------------------------------------
switch invcfg.regionPatt
    case 'prior'
        switch invcfg.type
            case 'weekly'
                patt = 'priorPatt_wi';
            case 'monthly'
                patt = 'priorPatt';
        end
        
    case 'flat'
        patt = 'flatPatt';
        
    otherwise
        error('Use "prior" or "flat" within-region pattern')
end


% Define Greens array
% ----------------------------------------------------------------------
% For 2011-13, resp should have dimensions 25x1096x3 (region x day x site).
% Greens will have dimensions 25x1096x3x36 (reg x day x site x month) for
% the monthly inversion, and 25x1096x3x144 (reg x day x site x week) for
% the weekly inversion.
% Greens may be diagonal in its 2nd dimension (see note below).
% NOTE, the month of source for a particular day is given by 
%       day2mon (each daily observation is influenced by sources 
%       from only one particular month - i.e. day2mon).
%       HOWEVER, for the weekly inversion, a daily observation may be
%       influenced by sources from a particular week or two consecutive
%       weeks (it is influenced by the past 4 days, which may be part of
%       two different weeks).
% UPDATE, the monthyl inversion now uses the same scheme as the weekly one
%         (2 months of influence are now permitted).
switch invcfg.type
    case 'monthly'
        greens13h = zeros(nreg,nday,nsites,nmon);
        greens15h = zeros(nreg,nday,nsites,nmon);
        for d = 1:nday
            % a few days early in a month are also influenced by the
            % previous month, consider this here.
            relD = relDays{d} - firstDay + 1;
            relMon = day2mon(relD(relD>0));
            
            % weight each month's contribution by its representation 
            % in relMon (ie. weight the corresponding equation)
            relMonU = unique(relMon);
            if length(relMonU) == 1
                % obs influenced by only one month
                greens13h(:,d,:,relMonU) = resp13h.(patt)(:,d,:);
                greens15h(:,d,:,relMonU) = resp15h.(patt)(:,d,:);
            else
                % obs influenced by two months
                for i = 1:length(relMonU)
                    weight = length(find(relMon==relMonU(i)))/length(relMon);
                    greens13h(:,d,:,relMonU(i)) = resp13h.(patt)(:,d,:)*weight;
                    greens15h(:,d,:,relMonU(i)) = resp15h.(patt)(:,d,:)*weight;
                end
            end
        end
        
    case 'weekly'
        greens13h = zeros(nreg,nday,nsites,nweek);
        greens15h = zeros(nreg,nday,nsites,nweek);
        for d = 1:nday
            % find the week(s) that influence this day
            % (on a weekly basis we can no longer neglect the fact that
            % about 50% of days are influenced by two weeks' sources)
            relD = relDays{d} - firstDay + 1;
            relWeek = day2week(relD(relD>0));
            
            % weight each week's contribution by its representation 
            % in relWeek (ie. weight the corresponding equation)
            relWeekU = unique(relWeek);
            if length(relWeekU) == 1
                % obs influenced by only one week
                greens13h(:,d,:,relWeekU) = resp13h.(patt)(:,d,:);
                greens15h(:,d,:,relWeekU) = resp15h.(patt)(:,d,:);
            else
                % obs influenced by two weeks
                for i = 1:length(relWeekU)
                    weight = length(find(relWeek==relWeekU(i)))/length(relWeek);
                    greens13h(:,d,:,relWeekU(i)) = resp13h.(patt)(:,d,:)*weight;
                    greens15h(:,d,:,relWeekU(i)) = resp15h.(patt)(:,d,:)*weight;
                end
            end
%             disp('TEST !!!')
%             greens13h(:,d,:,iweek(d)) = resp13h.(patt)(:,d,:);
%             greens15h(:,d,:,iweek(d)) = resp15h.(patt)(:,d,:);
        end
end


% Include 13h and/or 15h observations
% ----------------------------------------------------------------------
if obscfg.include13h && obscfg.include15h
    greens = cat(2,greens13h,greens15h);
elseif obscfg.include13h
    greens = greens13h;
elseif obscfg.include15h
    greens = greens15h;
else
    error('Choose 13h and/or 15h observations')
end
clear greens13h greens15h    
    
    
% % Permute dimensions to 731x3x25x24 (monthly 2011-2012 @15h inversion)
% greens = permute(greens,[2 3 1 4]);
% tmp = permute(greens,[1 2 4 3]);
% frech = reshape(tmp,ndat,nsrc);


% Build Frechet matrix
% ----------------------------------------------------------------------
% The Frechet matrix is just the reshaped Green's array, with
% dimensions ndat x nsrc ([731x3] x [25x24] for monthly 2011-12 @15h inv)
% NOTE, need to permute array due to subtleties of reshape function 
%       for higher-dimension arrays.
greens = permute(greens,[2 3 4 1]);
frech = reshape(greens,ndat,nsrc);

% Add the station offsets (ie. last 3 unknowns)
for s = 1:nsites
    I = (s-1)*nday*ntimes+1:s*nday*ntimes;
    frech(I,nsrc+s) = 1.0;
end

% Remove equations (=rows) with unavailable observations (obs = NaN)
frech = frech(iavail,:);


% Temporal smoothing of fluxes
% ----------------------------------------------------------------------
w2wdev_land = invcfg.smoothWland;
w2wdev_oce = invcfg.smoothWoce;

if (w2wdev_land>0 || w2wdev_oce>0) && strcmp(invcfg.type,'weekly')
    
    nweek = 48*(diff(invcfg.period)+1);
    ireg = 1:25; % regions to be smoothed
    nreg = length(ireg);
    lastrows = zeros(nreg*nweek,nsrc+nsites);
    
    for j = 1:nreg
        r = ireg(j);
        idx = ((r-1)*nweek+1):(r*nweek); % week index
        for w = 1:nweek-1
            lastrows(w+(j-1)*nweek,idx(w:w+1)) = [1 -1];
        end
        lastrows(nweek+(j-1)*nweek,idx([nweek 1])) = [1 -1];
    end
    
    frech = [frech; lastrows];
    
end


% Annual constraint for Open Ocean (OO) and Australia (Oz)
% ----------------------------------------------------------------------
annunc_pc = invcfg.OzOOAnn_pc; % per cent

if (annunc_pc > 0) && invcfg.oceanprior
    
    % number of sources per region and year
    switch invcfg.type
        case 'weekly'
            ns = 48;
        case 'monthly'
            ns = 12;
    end
    
    % number of years x Oz and OO regions ([16 20:25])
    ireg = [16 20:25];
    nreg = length(ireg);
    nyr = (diff(invcfg.period)+1);
    lastrows = zeros(nreg*nyr,nsrc+nsites);
    
    for j = 1:nreg
        r = ireg(j);
        for y = 1:nyr
            idx = ((r-1)*nyr*ns+(y-1)*ns+1):((r-1)*nyr*ns+y*ns); % source index
            lastrows((j-1)*nyr+y,idx) = ones(1,ns);
        end
    end
    
    frech = [frech; lastrows];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%