%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This routine calls the NZ inversion (regNZinv.m) routine, applying
% configurations according to the config file.
% 
% The reduced chi^2 value is computed after the inversion, so that 
% observational uncertainty can be adjusted in order to achieve chi^2=1.
% This is optional and configured in the config file.
% 
% Author: Kay Steinkamp
% Date: Jan/Feb 2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define local environment
% ----------------------------------------------------------------------
addpath '~/ResearchProjects/NZinv_GitHub'

% set absolute path to configuration file
configfile = '~/ResearchProjects/NZinv_GitHub/config';

fclose('all');


% set up a few global variables
% ----------------------------------------------------------------------
global monfrac monspan weekfrac weekspan monthnm weeknm 
global x2c ystr Nday dayrange

monthnm = {'jan','feb','mar','apr','may','jun',...
    'jul','aug','sep','oct','nov','dec'};
weeknm = cell(1,4*length(monthnm));
[weeknm{:}] = deal('w'); 
weeknm = cellstr([[weeknm{:}].',num2str((1:48)')])';
       
% lengths months and years 2011 to 2013
monlen11 = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
monlen12 = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
monlen13 = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
daysperyr11 = sum(monlen11);
daysperyr12 = sum(monlen12);
daysperyr13 = sum(monlen13);

% set each month to have 4 weeks (to ease later aggregation)
for m = 1:12
    W = 1+(m-1)*4:m*4;
    weeklen11(W) = monlen11(m)/4; %#ok<SAGROW>
    weeklen12(W) = monlen12(m)/4; %#ok<SAGROW>
    weeklen13(W) = monlen13(m)/4; %#ok<SAGROW>
end

% convert mole fraction to concentration using typical values for P, T
x2c = (44.01/8.31) * 97000 / 286;


% get the configuration parameters
% ----------------------------------------------------------------------
disp(' ')
disp('------------ Start of inversion procedure ----------------------')

[invcfg, obscfg, outcfg] = read_config(configfile);

% go to inversion directory
cd(invcfg.workdir); 

% finalize file names depending on inversion period
ystr = [num2str(invcfg.period(1)-2000) num2str(invcfg.period(2)-2000)];
invcfg.Fprior = [invcfg.Fprior '_' ystr '.mat'];
invcfg.FGreensFun = [invcfg.FGreensFun '_' ystr '.mat'];
obscfg.Fsynobs = [obscfg.Fsynobs '_' ystr '.mat'];
obscfg.Fobs = [obscfg.Fobs '_' ystr '.mat'];
obscfg.FrespFoss = [obscfg.FrespFoss '_' ystr '.mat'];
obscfg.FwtBsl = [obscfg.FwtBsl '_' ystr '.mat'];

% day range depending on inversion period
dayrange = {['01-Jan-20' ystr(1:2)]; ['31-Dec-20' ystr(3:4)]}; 

% put together weekfrac, weekspan as well as monfrac, monspan depending
% on inversion period
% *frac = fraction of year for weeks/months
% *span = span of days for weeks/months
years = str2double(ystr(1:2)):str2double(ystr(3:4));
monfrac = []; monspan = []; 
weekfrac = []; weekspan = [];
Nday = 0;
for y = years
    eval(['monfrac = [monfrac monlen' num2str(y) ...
        '/daysperyr' num2str(y) '];'])
    eval(['weekfrac = [weekfrac weeklen' num2str(y) ...
        '/daysperyr' num2str(y) '];'])   
    eval(['monspan = [monspan [cumsum(monlen' num2str(y) ...
        ')-monlen' num2str(y) '; cumsum(monlen' num2str(y) ')]+Nday];']);
    eval(['weekspan = [weekspan [cumsum(weeklen' num2str(y) ...
        ')-weeklen' num2str(y) '; cumsum(weeklen' num2str(y) ')]+Nday];']); 
    eval(['Nday = Nday + daysperyr' num2str(y) ';'])
end


% call inversion routine
% ----------------------------------------------------------------------
[data, sources] = regNZinv(invcfg, obscfg, outcfg);


% calculate reduced chi^2 value
% ----------------------------------------------------------------------
srcupdate = sources.postV(1:sources.nsrc)-sources.priorV(1:sources.nsrc);

ndat = data.neqns - data.nextra;   
mismatch = data.postCV(1:ndat) - data.priorCV(1:ndat);

sumsrc = sum(srcupdate.^2 ./ sources.prioruncV(1:sources.nsrc).^2);
sumdat = sum(mismatch.^2 ./ data.priorCuncV(1:ndat).^2);
% devide by DoF = Ndat+Nsrc-Nsrc
chi2 = (sumdat + sumsrc)/ndat;

fprintf('\nReduced chi^2: %4.2f\n',chi2);
str = sprintf('NZ %d-%d inversion complete',...
    invcfg.period(1),invcfg.period(2));
disp(' ');
disp(['----------- ' str ' --------------------']);
disp(' ');


% optimize chi^2 value
% ----------------------------------------------------------------------
if invcfg.optchi2
    step = int32(0);
    nstep = 0;
    obscfg.scale = obscfg.iniscale;
    while (chi2 < 0.99 || chi2 > 1.01) && (nstep < 10)
        nstep = nstep + 1;
        if chi2 < 0.99
            % chi2 too small, ie. data unc. should be decreased
            obscfg.scale = obscfg.scale * (1.0 - 0.8*(1.0 - chi2));
        else
            % mchi2 too large, ie. data unc. should be increased
            if chi2 >= 4.0
                error('doinv: chi2 too large (>= 4.0) !')
            elseif chi2 >= 2.0
                obscfg.scale = obscfg.scale * (1.0 + 0.3*(chi2 - 1.0));
            else
                obscfg.scale = obscfg.scale * (1.0 + (chi2 - 1.0));
            end
        end
        
        % run inversion again
        close all
        [data, sources] = regNZinv(invcfg, obscfg, outcfg);
        
        % calculate model mean chi^2 value as before
        srcupdate = sources.postV(1:sources.nsrc)-sources.priorV(1:sources.nsrc);
        mismatch = data.postCV(1:ndat) - data.priorCV(1:ndat);
        sumsrc = sum(srcupdate.^2 ./ sources.prioruncV(1:sources.nsrc).^2);
        sumdat = sum(mismatch.^2 ./ data.priorCuncV(1:ndat).^2);
        chi2 = (sumdat + sumsrc)/ndat; % devide by DoF = Ndat+Nsrc-Nsrc
        
        fprintf('\nReduced chi^2: %4.2f\n',chi2);
    end
end

data.chi2 = chi2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%