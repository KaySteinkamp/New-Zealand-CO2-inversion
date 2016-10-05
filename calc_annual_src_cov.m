function sources = calc_annual_src_cov(sources, invtype, toggle)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Computes the covariance matrix as well as source estimate for 
% group specified in the script (annual fluxes for all regions)
%
% The structure "sources" contains the source estimates, which are 
% supposed to be monthly fluxes in Tg CO2 yr-1.
%
% "nreg" is the number of surface regions.
% The string "toggle" switches between prior and posterior:
% toggle='prior' or toggle='post'
%
% Results are appended to the structure "sources"
%
% Author: Kay Steinkamp
% Date: Jan 2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get needed parameters
global monfrac ystr

nspatial = sources.nsrc;
nreg = sources.nreg;
nmon = length(monfrac);
if strcmp(invtype,'weekly')
    imon = reshape(repmat(sources.imonth,[1 nreg]),nreg*nmon,1);
else
    imon = sources.month(1:nspatial);
end
years = str2double(ystr(1:2)):str2double(ystr(3:4));
monperyr = nmon/length(years);


% switch between prior and posterior
% ----------------------------------------------------------------------
switch toggle
    case 'prior'
        if strcmp(invtype,'weekly')
            S = reshape(sources.priorM,nmon*nreg,1);
            C = diag(reshape(sources.prioruncM,nmon*nreg,1).^2);
        else
            S = sources.priorV(1:nspatial);
            C = diag(sources.prioruncV(1:nspatial).^2);
        end
    case 'post'
        if strcmp(invtype,'weekly')
            S = reshape(sources.postM,nmon*nreg,1);
            C = sources.postcovM;
        else
            S = sources.postV(1:nspatial);
            C = sources.postcov(1:nspatial,1:nspatial);
        end
    otherwise
        error('Specify either "prior" or "post"!')
end


% convert sources from Tg CO2 yr-1 to Tg CO2 month-1
% ----------------------------------------------------------------------
mfrac = monfrac(imon)';
mS = S .* mfrac;
mC = diag(mfrac) * C * diag(mfrac);


% calculate annual sources & covariances
% ----------------------------------------------------------------------
for r = 1:nreg
    
    i{r} = false(size(mS));
    i{r}((r-1)*nmon+1:r*nmon) = true;
    
    for y = 1:length(years)
        
        I = (i{r} & (imon>(y-1)*monperyr) & (imon<=y*monperyr));
    
        % annual mean sources for each year (Tg CO2 yr-1)
        aS{y}(r) = sum(mS(I));

        % full covariance of annual mean sources
        for s = 1:r
            J = (i{s} & (imon>(y-1)*12) & (imon<=y*12));
            aC{y}(r,s) = sum(sum(mC(I,J)));
            aC{y}(s,r) = aC{y}(r,s);
        end
        
        % uncertainty (Tg CO2 yr-1)
        uS{y}(r) = sqrt(aC{y}(r,r));
        
    end
    
end


% check positive definitness
% (sensitive to highest eigenvalue)
% ----------------------------------------------------------------------
tol = 1e-10;
for y = 1:length(years)
    if any(eig(aC{y}) < -tol*max(eig(aC{y}))) 
        disp(['WARNING: Annual ',toggle,' covariance matrix not positive definit!']);
    end
end
    

% append results to sources structure
% ----------------------------------------------------------------------
sources.years = cellstr(num2str(years' + 2000))';
sources.([toggle 'A']) = aS;
sources.([toggle 'uncA']) = uS;
if strcmp(toggle,'post')
    sources.([toggle 'covA']) = aC;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%