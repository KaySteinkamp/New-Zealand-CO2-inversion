function sources = read_priorsources(invcfg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Reads monthly or weekly prior source estimates for the NZ 
% atmospheric CO2 inversion.
%
% In total, 25*24 + 1 = 601 sources are considered.
% -> The first (25*24) represent the monthly regional sources, which are 
%    read from a mat file and are based on the Biome-BGC (land) and 
%    ROMS-PISCES/NEMO (ocean) models. 
% -> The offset value (i.e. the last extra tracer) is predetermined in 
%    the config file, and should be set to zero with large uncertainty.
%
% Returned quantities are stored in the structure "sources".
% They include the names of sources, their index, values and variance.
% 
% Author: Kay Steinkamp
% Date: Jan 2014, updated in Mar 2014 to include weekly priors
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global x2c ystr weekfrac monfrac

% Read prior sources (in Tg CO2 yr-1) for all regions and months/weeks
% ----------------------------------------------------------------------
load([invcfg.workdir invcfg.Fprior], 'regprior');

nreg = length(regprior.regname);
imon = regprior.(['mon' ystr]);
nmon = length(imon);
if ~isequal(imon,(1:nmon)')
    error('Unexpected monthly indices')
end
month = repmat(imon,[1 nreg]);

switch invcfg.type
    case 'monthly'
        if strcmp(invcfg.ocpriorType,'PISCES')
            pest = regprior.scaleFac';
            punc = regprior.scaleFacUnc';
        elseif strcmp(invcfg.ocpriorType,'Takahashi')
            pest = regprior.scaleFac_Tak';
            punc = regprior.scaleFacUnc_Tak';
        else
            error('Invalid type of ocean prior: Takahashi or PISCES')
        end
        % Reshape reg*month matrix to vector
        month = reshape(month,nmon*nreg,1);
        week = []; week2mon = []; nweek = []; weekmidINmon = [];
        nt = nmon;
    case 'weekly'
        if strcmp(invcfg.ocpriorType,'PISCES')
            pest = regprior.scaleFac_wi';
            punc = regprior.scaleFacUnc_wi';
        elseif strcmp(invcfg.ocpriorType,'Takahashi')
            pest = regprior.scaleFac_Tak_wi';
            punc = regprior.scaleFacUnc_Tak_wi';
        else
            error('Invalid type of ocean prior: Takahashi or PISCES')
        end
        % Reshape reg*month and reg*week matrices to vectors
        % NOTE, each month is subdivided into 4 weeks
        weekmidINmon = regprior.weekmidINmon;
        nweek = length(weekmidINmon);
        for m = 1:nmon
            W = (1+(m-1)*4:m*4)';
            week(W,:) = repmat(W,[1 nreg]);
            week2mon(W) = m;
        end 
        month = reshape(month(week2mon,:),nweek*nreg,1);
        week = reshape(week,nweek*nreg,1);
        nt = nweek;
    otherwise
        error('Use "monthly" or "weekly" as inversion type')
end


% Switch priors on/off
if ~invcfg.landprior
    punc(:,1:16) = 1e8;
end
if ~invcfg.oceanprior
    punc(:,17:25) = 1e8;
end
if ~invcfg.OzOOprior
    punc(:,[16 20:25]) = 1e8;
end

nsrc = numel(pest);

% Reshape reg*month (or reg*week) matrix to vector
pestV = reshape(pest,nsrc,1);
puncV = reshape(punc,nsrc,1);

if strcmp(regprior.scaleFacUnit,'Multiples of 1 Tg CO2 yr-1')
    units = 'Tg CO2 yr-1';
end


% Offset concentration (last extra tracer(s)).
% Use some average T and P values to convert from ppm to g CO2 m-3
% ----------------------------------------------------------------------
offnm = {'Offset BHD';'Offset LAU';'Offset RBM'};
if length(offnm)~=length(invcfg.offset)
    error('Screwed up station offsets')
end
offestX = invcfg.offset;
offuncX = invcfg.offsetUnc;
offestC = offestX*1e-6 * x2c;
offuncC = offuncX*1e-6 * x2c;


% Store in sources structure
% ----------------------------------------------------------------------
sources.nreg = nreg;
sources.regarea = regprior.regarea;
sources.imonth = imon;
sources.nmon = nmon;
sources.week2mon = week2mon';
sources.nweek = nweek;
sources.weekmidINmon = weekmidINmon';
sources.weekfrac = weekfrac;
sources.monfrac = monfrac;
sources.name = [regprior.regname; offnm];
sources.nsrc = nsrc;
sources.nvar = nsrc + length(offestX);
sources.nt = nt;
sources.month = [month; zeros(length(offestX),1)];
sources.week = [week; zeros(length(offestX),1)];
sources.prior = pest;
sources.priorunc = punc;
sources.priorV = [pestV; offestC];
sources.prioruncV = [puncV; offuncC];
sources.unit = units;
sources.offestX = offestX;
sources.offuncX = offuncX;
sources.offestC = offestC;
sources.offuncC = offuncC;

                                  
% Check symmetry & positive definitness of prior covariance matrix
% ----------------------------------------------------------------------
% NOTE, this test is currently obsolete (prior cov matrix is diagonal)
% tol = 1e-9;
% pcov = sources.priorcov;
% for n = 1:size(pcov,3)
%     if any(any(abs(pcov(:,:,n) - pcov(:,:,n)') > tol))
%         error('make_priorsource: prior covariance matrix not symmetric!');
%     end
%     if any(eig(pcov(:,:,n)) < 0)
%         error('make_priorsource: prior covariance matrix not positive definit!');
%     end
% end


% calculate/plot correlation matrix
% ----------------------------------------------------------------------
% for n = 1:size(pcov,3)
%     pcorr = diag(diag(pcov(:,:,n)).^-0.5) * pcov(:,:,n) * ...
%             diag(diag(pcov(:,:,n)).^-0.5);
%     figure('renderer','zbuffer')
%     pcolor(pcorr); axis ij; shading flat
%     colorbar('peer',gca);
%     title('Prior Source Correlation Matrix')
%     %close
% end
% clear pcov


% Display some information
% ----------------------------------------------------------------------
if invcfg.chatty                                         
    fprintf('\nPrior sources have been read for %g-%g\n',...
        invcfg.period(1),invcfg.period(2));
    fprintf('Number of source regions  : %g\n',sources.nreg);
    fprintf('Number of months in period: %g\n',sources.nmon);
    if ~isempty(sources.nweek)
        fprintf('Number of weeks in period : %g\n',sources.nweek);
    end
    fprintf('Total number of unknowns  : %g\n',sources.nvar);
    fprintf('Range of prior estimates  : %5.2f to %5.2f Tg CO2 yr-1\n',...
        min(pestV),max(pestV));
    fprintf('Range of prior uncertainty: %5.2f to %5.2f Tg CO2 yr-1\n',...
        min(puncV),max(puncV));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%