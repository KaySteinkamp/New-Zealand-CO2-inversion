function [data, sources] = regNZinv(invcfg, obscfg, outcfg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Performs weekly/monthly Bayesian inversion of 2xdaily CO2 records.
% 
% Based on a set of regional sources with predetermined within-region
% spatial pattern.
%
% The input data is the observational response for each source (from 
% Lagrangian disperion backward modelling), a series of prior estimates 
% for each source (along with variances) and a record of observations 
% (daily afternoon CO2 measurements at BHD, LAU and RBM through 2011-13)
%
% There are no additional constraints on sources (or groups of sources).
%
% Along with the sources, a constant observational offset is estimated,
% which can be freely (prior uncertainty large) set by the inversion.
%
% The Green's functions, ie. observational response, are precalculated 
% by backward Lagrangian dispersion in the NAME model, using hourly
% meteorology from the NZLAM-12 model.
%
% The computational core of the code is the routine bayesinv.m which is
% a general-purpose Bayesian inverse solver.
%
% 
% Author: Kay Steinkamp (kay.steinkamp@niwa.co.nz)
% Date: Jan 2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global x2c

% Get prior source
% ----------------------------------------------------------------------
sources = read_priorsources(invcfg);


% Get observational data
% (sources needed as input for extra constraints)
% ----------------------------------------------------------------------
data = read_data(invcfg, obscfg, sources);

          
% Get model responses, i.e. the Green's functions, and
% build Frechet/Jacobi matrix
% ----------------------------------------------------------------------
frech = makefrech(invcfg, obscfg, data.nsites, data.ntimes, ...
                  data.ndat, data.iavail, sources);

if invcfg.chatty                                         
    disp(' ')
    disp('Frechet Matrix built')
end


% Compute temporal aggregates for prior sources 
% (annual/monthly)
% ----------------------------------------------------------------------
if strcmp(invcfg.type,'weekly')
    sources = calc_monthly_src_cov(sources, 'prior');
    if invcfg.chatty
        disp(' ')
        disp('Monthly prior source and covariance computed')
    end
end
sources = calc_annual_src_cov(sources, invcfg.type, 'prior');
if invcfg.chatty
    disp('Annual prior source and covariance computed')
end


% Compute "prior pseudo data" (multiply frechet matrix with prior sources)
% ----------------------------------------------------------------------
data.pseudoCV = frech * sources.priorV;
data.pseudoXV = (8.31/44) * 1e6 * (data.pseudoCV .* data.TV ./ data.PV);
data.pseudoX = nan(size(data.priorX));
data.pseudoX(data.iavail) = data.pseudoXV(1:length(data.iavail));


% Output prior source & data information
% ---------------------------------------------------------------------- 
outcfg = write_output(invcfg, obscfg, outcfg, sources, data, 'prior'); 

if invcfg.chatty                                         
    disp(' ')
    disp('Prior output written')
end

                     
% Call the inversion core (SVD based computation)
% ---------------------------------------------------------------------- 
if invcfg.chatty
    disp(' ')
    disp('Calling inversion core...')
    tic
end
% full posterior data covariance not needed, only variance,
% so set toggle (last argument of bayesinv) to 0.
[sources.postV, sources.postcov, ...
   data.postCV, data.postCuncV] = ...
                bayesinv(sources.priorV, diag(sources.prioruncV.^2), ...
                         data.priorCV, data.priorCuncV.^2, frech, 0);

data.postCuncV = sqrt(data.postCuncV);

if invcfg.chatty
    t = toc;
    disp(['Inversion core completed in ' num2str(round(t)) ' seconds'])
end


% Calculate CO2 mole fraction from concentration;
% also reshape to convenient form
% ----------------------------------------------------------------------
data.postXV = (8.31/44) * 1e6 * (data.postCV .* data.TV ./ data.PV);
data.postXuncV = (8.31/44) * 1e6 * (data.postCuncV .* data.TV ./ data.PV);

% Also reshape to convenient form
data.postC = nan(size(data.priorC));
data.postCunc = nan(size(data.priorCunc));
data.postC(data.iavail) = data.postCV(1:length(data.iavail));
data.postCunc(data.iavail) = data.postCuncV(1:length(data.iavail));

data.postX = nan(size(data.priorX));
data.postXunc = nan(size(data.priorXunc));
data.postX(data.iavail) = data.postXV(1:length(data.iavail));
data.postXunc(data.iavail) = data.postXuncV(1:length(data.iavail));


% Reshape source array to more convenient, non-vector, form
% ----------------------------------------------------------------------
sources.postuncV = sqrt(diag(sources.postcov));
sources.post = reshape(sources.postV(1:sources.nsrc), ...
                        sources.nt,sources.nreg);
sources.postunc = reshape(sources.postuncV(1:sources.nsrc), ...
                           sources.nt,sources.nreg);
                       
% convert offset back from concentration to mole fraction
sources.offestC(:,2) = sources.postV(sources.nsrc+1:sources.nvar);
sources.offuncC(:,2) = sources.postuncV(sources.nsrc+1:sources.nvar);
sources.offestX(:,2) = sources.offestC(:,2)/x2c*1e6;
sources.offuncX(:,2) = sources.offuncC(:,2)/x2c*1e6;


% Compute temporal aggregates for prior sources 
% (annual/monthly)
% ----------------------------------------------------------------------
if strcmp(invcfg.type,'weekly')
    sources = calc_monthly_src_cov(sources, 'post');
    if invcfg.chatty
        disp(' ')
        disp('Monthly posterior source and covariance computed')
    end
end
sources = calc_annual_src_cov(sources, invcfg.type, 'post');
if invcfg.chatty
    disp('Annual prior source and covariance computed')
end


% Output posterior source & data information
% ---------------------------------------------------------------------- 
write_output(invcfg, obscfg, outcfg, sources, data, 'post'); 

if invcfg.chatty                                         
    disp(' ')
    disp('Posterior output written')
end

         
% Do some quick plots of the inversion results and store them
% ----------------------------------------------------------------------
if outcfg.doplots
    makeplots(invcfg, outcfg, data, sources)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%