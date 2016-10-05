function [invcfg, obscfg, outcfg] = read_config(cfgfile)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Reads the control parameters for the NZ inversion.
%
% Parameters are read from a configuration file "cfgfile" and 
% returned in structures. 
%
% 
% Author: Kay Steinkamp
% Date: Jan 2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% **********************************************************************
% open the configuration file and read control parameters
% **********************************************************************
ison = @(str) strcmp(str,'.true.') || strcmp(str,'t');  %#ok<NASGU>

f10 = fopen(cfgfile,'r');


% 1) inversion configuration
% -----------------------------
dum = textscan(f10,'%d',1,'HeaderLines',10);
count = dum{1};

dum = textscan(f10,'%d %s = %s',count);
type = dum{1};
name = dum{2};
value = dum{3};

for i = 1:count
    switch type(i)
        case 0
            eval(['tmp.' name{i} '= ison(value{i});']);
        case 1
            eval(['tmp.' name{i} '= int32(str2double(value{i}));']);
        case 2
            eval(['tmp.' name{i} '= str2double(value{i});']);
        case 3
            eval(['tmp.' name{i} '=' value{i} ';']);        
        otherwise
            error('read_config: invalid data type in configuration file!')
    end
end

invcfg = tmp;
clear tmp;


% 2) handling of observations
% -----------------------------
dum = textscan(f10,'%d',1,'HeaderLines',3);
count = dum{1};

dum = textscan(f10,'%d %s = %s',count);
type = dum{1};
name = dum{2};
value = dum{3};

for i = 1:count
    switch type(i)
        case 0
            eval(['tmp.' name{i} '= ison(value{i});']);
        case 2
            eval(['tmp.' name{i} '= str2double(value{i});']);
        case 3
            eval(['tmp.' name{i} '=' value{i} ';']);
        otherwise
            error('read_config: invalid data type in configuration file!')
    end
end

obscfg = tmp;
clear tmp;


% 3) output configuration
% -----------------------------
dum = textscan(f10,'%d',1,'HeaderLines',3);
count = dum{1};

dum = textscan(f10,'%d %s = %s',count);
type = dum{1};
name = dum{2};
value = dum{3};

for i = 1:count
    switch type(i)
        case 0
            eval(['tmp.' name{i} '= ison(value{i});']);
        case 1
            eval(['tmp.' name{i} '= int32(str2double(value{i}));']);
        case 2
            eval(['tmp.' name{i} '= str2double(value{i});']);
        case 3
            eval(['tmp.' name{i} '=' value{i} ';']);
        otherwise
            error('read_config: invalid data type in configuration file!')
    end
end

outcfg = tmp;
clear tmp;

fclose(f10);


% **********************************************************************
% make sure that all needed parameters are actually present
% **********************************************************************
% needed parameters for inversion
invneeded = {...
'workdir', ...        % working directory
'type', ...           % type of inversion (monthly or weekly)
'Fprior', ...         % filename of prior sources 
'Ffemi', ...          % filename of fossil fuel emissions 
'landprior', ...      % flag to include land prior
'oceanprior', ...     % flag to include ocean prior
'ocpriorType', ...    % choose 'Takahashi' or 'PISCES'
'OzOOprior', ...      % flag to include Australian and OpenOcean prior
'OzOOAnn_pc', ...     % percentage of Australian and OpenOcean annual unc
'regionPatt', ...     % type of within-region flux pattern
'smoothWland', ...    % Gaussian temporal smoother for weekly land fluxes
'smoothWoce', ...     % Gaussian temporal smoother for weekly ocean fluxes
'offset', ...         % prior value of constant offset (ppm)
'offsetUnc', ...      % prior uncertainty of constant offset (ppm)
'period', ...         % period of inversion
'FGreensFun', ...     % filename of Green's Function for 2011-12
'optchi2', ...        % flag to adjust data uncertainty to achieve chi2=1
'chatty'};            % be chatty about progress

% needed parameters for observations
obsneeded = {...
'syntheticObs', ...   % flag to use synthetic observations
'synthDiCy', ...      % flag to use synthetic observations with diurnal cycle from BiomeBGC
'presubFoss', ...     % flag to presubtract response to foss. em. from obs.
'includeBHD', ...     % flag to include observations from BHD
'includeLAU', ...     % flag to include observations from LAU
'includeRBM', ...     % flag to include observations from RBM
'Fsynobs', ...        % filename of synthetic observations 
'include13h', ...     % flag to include 13-14h LT observations
'include15h', ...     % flag to include 15-16h LT observations
'FrespFoss', ...      % filename of obs. response to fossil emissions 
'synobsUnc', ...      % uncertainty for synthetic observations 
'Fobs', ...           % filename of observations 
'useWtBsl', ...       % flag to use weighted baseline (BHD/TF5 for S/N) 
'BslBiasSig', ...     % applied baseline bias, in multiples of baseline st-dev 
'adjustCCL', ...      % flag to adjust data to NOAA CCL using Aniwaniwa suite (Oct 2014) 
'smoothObs', ...      % flag to smooth obs to match up source/obs timescale 
'FwtBsl', ...         % filename of weighted baseline 
'minUnc', ...         % minimum data uncertainty, in ppm 
'uniNoiseX', ...      % uniform noise, in ppm 
'biasX', ...          % observation bias, in ppm 
'iniscale'};          % initial scaling factor of observational uncertainty 

% needed parameters for output
outneeded = {...
'writedir', ...       % path to output directory    
'dosources', ...      % flag to write out sources (prior & posterior)
'dodata', ...         % flag to write out posterior data
'doplots'};           % flag to make plots (and store them)

% now check
if ~(all(isfield(invcfg,invneeded)) && all(isfield(obscfg,obsneeded)) ...
        && all(isfield(outcfg,outneeded)))
    error('read_config: required parameter(s) not present in config file!')
end

% **********************************************************************
% copy the employed configuration to output folder
% **********************************************************************
if ~exist([invcfg.workdir outcfg.writedir], 'dir')
    mkdir([invcfg.workdir outcfg.writedir])
end
copyfile(cfgfile,[invcfg.workdir outcfg.writedir 'configuration'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%