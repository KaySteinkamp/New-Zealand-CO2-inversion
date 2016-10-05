function [src2, csrc2, dat2, cdat2] = ...
                           bayesinv(src, csrc, dat, cdat, pdm, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script is a general-purpose Bayesian inverse solver.
%
% Bayesian inversion is performed via SVD
% 
% It is based on the IDL wrapper for transcom 3, Level 1 inversions, 
% which originally was written by:
% Peter Rayner, Kevin Gurney, Rachel Law, October 2000
%
% The code was extended in 2008 by Kay Steinkamp in order to:
% (1) cope with under-determined systems (i.e. more sources than data)
% (2) allow for non-diagonal covariance matrices (both data and sources)
%
% ----------------------------------------------------------------------
% In more detail:
% this routine solves the Bayesian least squares problem
% pdm * src = dat 
% subject to prior estimates src and dat, with covariances cdat and csrc
%
% It does this by minimizing the quantity
% |dat2 - dat|^2 + |src2 - src|^2 
% weighted by their inverse covariances
% 
% The equations are taken from Tarantola and Valette (1982) in their
% linear forms. The equations are
%
% src2 = (pdm~ cdat^-1 pdm + csrc^-1)^-1*
%        (pdm~ cdat^-1 dat + csrc^-1 src)   (1)
%
% dat2 = pdm * src2                         (2)
%
% csrc2 = (pdm~ cdat^-1 pdm + csrc^-1)^-1   (3)
% 
% cdat2 = pdm csrc2 pdm~                    (4)
%
% pdm is "partial derivative matrix", i.e. Frechet (Jacobian) matrix
% src is model parameter vector, i.e. the sources
% dat is data vector
% initial "c" means covariance,
% suffix "2" means posterior value,
% no suffix means prior value and
% "~" indicates transpose
% 
% ----------------------------------------------------------------------
% The routine solves the Bayesian inversion problem (as depicted above)
% via an SVD solution of an augmented nonBayesian problem.
% 
% The SVD is performed on the studentized pdm: 
% pdm_st = cdat^-0.5 pdm csrc^0.5
% U*S*V~ = pdm_st
%
% The above equations can then be replaced by:
%
% src2 = (csrc^0.5 V S(S^2+1)^-1 U~ cdat^-0.5) dat +
%        (1 - csrc^0.5 V S^2(S^2+1)^-1 V~ csrc^-0.5) src   (1*)
%
% incr = cdat^0.5 (U S^2(S^2+1)^-1 U~ - 1) cdat^-0.5      (auxiliary)
% dat2 = (incr + 1) dat -
%        incr cdat^0.5 U S V~ csrc^-0.5 src                (2*)
%
% csrc2 = csrc^0.5 (1 - V S^2 (S^2+1)^-1 V~) csrc^0.5      (3*)
%
% cdat2 implemented, but not verified yet (never needed)
% 
% ----------------------------------------------------------------------
% In case of diagonal covariance matrices (cdat | csrc), the code is
% significantly faster than with "full" matrices.
% cdat can alternatively provided as vector, this routine then assumes
% that the vector contains the data variances.
%
% ----------------------------------------------------------------------
% Example for dimensions of prior arrays (row major notation):
% src:  (27x1)              
% csrc: (27x27) diagonal/non-diagonal
% dat:  (77x1)
% cdat: (77x77) diagonal/non-diagonal
% pdm:  (77x27)
%
% Example for dimensions of posterior arrays:
% src2: (27x1)
% csrc2:(27x27) non-diagonal
% dat2: (77x1)
%
% ----------------------------------------------------------------------
% Author: Kay Steinkamp
% Date: Nov 2008
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set optional parameters
doplot = false;   % default value
docdat2 = true;   % default value
if nargin > 5
    docdat2 = varargin{1};
    if nargin > 6
        doplot = varargin{2};
    end
end

% flag whether cdat is diagonal, 
% or if given as vector of variances.
% transform it into a sparse matrix in both cases
cdatvec = false;  % default
if xor(size(cdat,1)==1, size(cdat,2)==1)
    % cdat is given as vector
    cdatvec = true;
elseif isempty(find(cdat - diag(diag(cdat)), 1))
    % cdat is given as diagonal matrix.
    % transform it into a vector
    cdat = diag(cdat);
    cdatvec = true;
end

% flag whether csrc is diagonal
if isempty(find(csrc - diag(diag(csrc)), 1))
    csrcdiag = true;
else
    csrcdiag = false;
end

% calculate residual (in terms of prior)
resid = dat - pdm*src;

% studentize residuals: devide by data uncertainty
if cdatvec
    resid = resid./cdat.^0.5;
else
    resid = sqrtm(inv(cdat)) * resid;
end

% studentize pdm: scale columns by (data uncertainty)^-1 
%                 and rows by prior uncertainty
if cdatvec && csrcdiag
    pdm = sparse(1:length(cdat),1:length(cdat),1./cdat.^0.5) * pdm * csrc.^0.5;
elseif csrcdiag
    pdm = sqrtm(inv(cdat)) * pdm * csrc.^0.5;
elseif cdatvec
    pdm = sparse(1:length(cdat),1:length(cdat),1./cdat.^0.5) * pdm * sqrtm(csrc);
end

% perform singular value decomposition
% NOTE: roles of rpdm/lpdm are exchanged compared to original IDL code
% [lpdm,spdm,rpdm] = svd(pdm,0);
[lpdm,spdm,rpdm] = svd(pdm,'econ');

% calculate the projections onto the singular vectors
picard = lpdm' * resid;

% check Discrete Picard Condition (DPC)
% ie. abs(U'*dat_weighted) should decay faster (on average) 
% to zero than the singular values (ie. diag(S))
if doplot
    sing = diag(spdm);
    xcoord = 1:size(spdm,1);
    plot(xcoord,sing,'r+');
    hold on
    plot(xcoord,abs(picard),'g*');
    hold off
    ylim([0 20])
end

% scale the projections by the Bayesian singular value
sinv = diag(spdm)./(diag(spdm).^2 + 1);
% NOTE: sinv is somewhat similar to the inverse of spdm, but is related 
%       to the nonBayesian-Bayesian transition (studentize)
picard = picard .* sinv;

% use the right s-vecs to construct the source increment
src2 = rpdm * picard;
% scale by std dev to get unstudentized residual
if csrcdiag
    src2 = (csrc.^0.5)'*src2 + src;
else
    src2 = sqrtm(csrc)'*src2 + src;
end

% now use singular value property to reconstruct data
picard = picard .* diag(spdm);
dat2 = lpdm*picard - resid; clear lpdm
% and unstudentize
if cdatvec
    % ssdat is sparse matrix with standard devs on main diagonal
    ssdat = sparse(1:length(cdat),1:length(cdat),cdat.^0.5);
    dat2 = ssdat*dat2 + dat;
else
    dat2 = sqrtm(cdat)*dat2 + dat;
end

% now posterior covariance
% fix nonbayesian to bayesian singular values
sfix = diag(diag(spdm)./(diag(spdm).^2 + 1).^0.5);
rpdm = rpdm * sfix;
% compute reduction
red = eye(size(csrc)) - rpdm*rpdm';

% source covariance
if csrcdiag
    csrc2 = (csrc.^0.5)' * red * csrc.^0.5;
else
    csrc2 = sqrtm(csrc)' * red * sqrtm(csrc);
end

% data posterior covariance
% Note that the following has not been verified yet
if docdat2
    if cdatvec
        cdat2 = (ssdat * pdm) * red * (ssdat * pdm)';
    else
        cdat2 = (sqrtm(cdat) * pdm) * red * (sqrtm(cdat) * pdm)';
    end
else
    % return only main diagonal (i.e. variances) 
    % if full covariance is not desired.
    % the goal is to compute the main diagonal of A*B 
    % without building the full matrix (A*B) first.
    if cdatvec
        A = ssdat * pdm;
    else
        A = sqrtm(cdat) * pdm;
    end
    clear pdm ssdat  % free up some memory
    B = red * A';
    cdat2 = zeros(size(A,1),1);
    for i = 1:size(A,1)
        cdat2(i) = A(i,:) * B(:,i); % result is a scalar
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%