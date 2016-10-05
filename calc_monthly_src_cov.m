function sources = calc_monthly_src_cov(sources, toggle)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Computes the covariance matrix as well as source estimate for 
% group specified in the script (monthly fluxes for all regions)
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
% Date: Mar 2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get needed parameters
global weekfrac monfrac

nspatial = sources.nsrc;
nreg = sources.nreg;
imon = sources.month(1:nspatial);
iweek = sources.week(1:nspatial);
nmon = sources.nmon;
nweek = sources.nweek;
week2mon = sources.week2mon;
if ~isequal(imon,week2mon(iweek))
    error('Temporal Information screwed up')
end


% switch between prior and posterior
% ----------------------------------------------------------------------
switch toggle
    case 'prior'
        S = sources.priorV(1:nspatial);
        C = diag(sources.prioruncV(1:nspatial).^2);
    case 'post'
        S = sources.postV(1:nspatial);
        C = sources.postcov(1:nspatial,1:nspatial);
    otherwise
        error('Specify either "prior" or "post"!')
end


% convert sources from Tg CO2 yr-1 to Tg CO2 week-1
% ----------------------------------------------------------------------
wS = S .* weekfrac(iweek)';
wC = diag(weekfrac(iweek)) * C * diag(weekfrac(iweek));


% calculate monthly mean sources & covariances
% ----------------------------------------------------------------------
for r = 1:nreg
    
    i{r} = false(size(wS));
    i{r}((r-1)*nweek+1:r*nweek) = true;
    
    for m = 1:nmon
        
        I = ( i{r} & (imon == m) );
    
        % monthly mean sources (Tg CO2 month-1)
        mS(m,r) = sum(wS(I));

        % full covariance for all reg/mon pairs
        for s = 1:r
            for n = 1:m
                J = ( i{s} & (imon == n) );
                mC(m,r,n,s) = sum(sum(wC(I,J)));
                mC(n,s,m,r) = mC(m,r,n,s);
            end
        end
        
        % uncertainty (Tg CO2 month-1)
        uS(m,r) = sqrt(mC(m,r,m,r));
        
    end
    
end


% reshape covariance 4D array into matrix
% ----------------------------------------------------------------------
mC = reshape(mC,nreg*nmon,nreg*nmon);


% convert sources from Tg CO2 month-1 to Tg CO2 yr-1
% ----------------------------------------------------------------------
mfrac = repmat(1./monfrac(sources.imonth)',[1 nreg]);
mmfrac = diag(reshape(mfrac,nreg*nmon,1));
mS = mS .* mfrac;
uS = uS .* mfrac;
mC = mmfrac * mC * mmfrac;


% check positive definitness
% (sensitive to highest eigenvalue)
% ----------------------------------------------------------------------
tol = 1e-10;
tmp = eig(mC);
if any(tmp < -tol*max(tmp))
    disp(['WARNING: monthly ',toggle,...
        ' covariance matrix may not be positive definit!']);
    disp(['Smallest eigenvalue ' num2str(min(tmp))])
%     if any(tmp < -tol*1e5*max(tmp))
%         error(['Monthly ',toggle,...
%             ' covariance matrix not positive definit!']);
%     end
end
    

% append results to sources structure
% ----------------------------------------------------------------------
sources.([toggle 'M']) = mS;
sources.([toggle 'uncM']) = uS;
if strcmp(toggle,'post')
    sources.([toggle 'covM']) = mC;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%