function [tsmon,tspresub] = read_responses(greenfile,nreg,elapm, ...
                                           npresub,data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script reads the responses from netCDF file for all available
% sites, i.e. those which were included in the forward model runs.
% 
% Then it picks those responses belonging to the sites used in the
% inversion (kept in "respindex").
%
% It returns the responses for all tracers (CO2 and presub) over the
% model period (elapm months), for all months pulsed, and for all
% regions (land and ocean). 
%
% For the script to work several netCDF routines have to be accessable.
%
% Author: Kay Steinkamp
% Date: Nov 2008
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get needed parameters
global monperyr

% **********************************************************************
% open and read greens function netCDF file
% **********************************************************************
ncid = netcdf.open(greenfile, 'NC_NOWRITE');

[~,nvars,~,~] = netcdf.inq(ncid);

% assign netCDF variables
for i = 0:nvars-1
    [varname,~,~,~] = netcdf.inqVar(ncid,i);
    switch varname
        case 'tracer_data_month'
            % response to CO2 tracers
            ncmon = permute(netcdf.getVar(ncid,i,'double'),[4 3 2 1]);
        case 'presub_data_month'
            % response to presub tracers
            ncpre = permute(netcdf.getVar(ncid,i,'double'),[3 2 1]);
        otherwise
            disp(['Skipping netcdf variable: ' varname])
    end
end
netcdf.close(ncid)
clear ncid varname nvars 

% nc = netcdf(greenfile, 'nowrite'); 
% 
% vars = var(nc);   
% 
% % assign netCDF variables
% for i = 1:length(vars)
%     switch name(vars{i})
%         case 'tracer_data_month'
%             ncmon = vars{i};    % response to CO2 tracers
%         case 'presub_data_month'
%             ncpre = vars{i};    % response to presub tracers
%         otherwise
%             error('read_responses.m: invalid name of netCDF variable')      
%     end
% end
% disp(greenfile)
%
% % close data netCDF file
% close(nc);

% squeeze out empty dimensions
tspresub = squeeze(ncpre);
tsmon = squeeze(ncmon);

clear ncpre ncmon

% check order of dimensions
nrespstat = size(tsmon,4);
% nrespstat = 245;     % number of stations available in responses

if ~isequal(size(tsmon),[elapm monperyr nreg nrespstat])
    error('read_responses: tracer response data in wrong order')
elseif ~isequal(size(tspresub),[elapm npresub nrespstat])
    error('read_responses: presub response data in wrong order')
end

% **********************************************************************
% select the responses we want (those associated to list of obs. sites)
% **********************************************************************
% Note: "respindex" translates the number of each site into its 
%       corresponding response function index (length(respindex)==nobs)
% Note: due to Matlab's indexing scheme, the following 2 lines are not
%       doing the right thing (bug in Matlab?). Thus a workaround is
%       applied involving permutation of the array dimensions
% tsmon = tsmon(:,:,:,respindex);
% tspresub = tspresub(:,:,respindex);
%
% Workaround:
% dum = permute(tsmon(:,:,:,:),[4 1 2 3]);
% tsmon = ipermute(dum(respindex,:,:,:),[4 1 2 3]);
% dum = permute(tspresub(:,:,:),[3 1 2]);
% tspresub = ipermute(dum(respindex,:,:),[3 1 2]);
%
% Update: in the Matlab 2011b version the bug doesn't exist anymore.
% **********************************************************************
% for the LSCE-LMDZ model the list of stations is different, so we need
% to assign the responses to the station names
% **********************************************************************
if isempty(strfind(greenfile,'LSCE'))
    tsmon = tsmon(:,:,:,data.respindex);
    tspresub = tspresub(:,:,data.respindex);
else
    % read LSCE site names and choose correct response indices
    names = lower(cellstr(ncread(greenfile,'site_name')'));
    LSCErespindex = zeros(data.nsites,1);
    for i = 1:data.nsites
        dname = strtok(data.sitename{i},'_');
        j = find(strcmp(dname,names));
        if ~isempty(j)
            LSCErespindex(i) = j;
        else
            % consider deviations in name scheme
            switch dname
                case 'bal'
                    j = find(strcmp('bal1',names));
                case 'car040'
                    j = find(strcmp('car_04000',names));
                case 'car050'
                    % car_05000 not in dataset, use car_04000 instead
                    j = find(strcmp('car_04000',names));                    
                case 'car060'
                    j = find(strcmp('car_06000',names));   
                case 'frd040'
                    j = find(strcmp('frd',names));            
                case 'ljo'
                    j = find(strcmp('ljo_w',names));            
                case 'gsn'
                    % Korean station gsn not in dataset, use tap instead
                    % (also in Korea, same longitude, 3° further north)
                    j = find(strcmp('tap',names));        
                case 'mhdcbc'
                    % Choose either mhd or mhd_w
                    j = find(strcmp('mhd_w',names));            
                case 'mhdrbc'
                    % Choose either mhd or mhd_w
                    j = find(strcmp('mhd_w',names));
                case 'stmebc'
                    j = find(strcmp('stm',names));
                case 'wpo000'
                    j = find(strcmp('wpo_00',names));    
                case 'wpon05'
                    % Western Pacific cruise not completely included,
                    % choose neighboring station instead.
                    j = find(strcmp('wpo_n10',names));
                case 'wpon10'
                    j = find(strcmp('wpo_n10',names));
                case 'wpon15'
                    % Western Pacific cruise not completely included,
                    % choose neighboring station instead.
                    j = find(strcmp('wpo_n20',names));
                case 'wpon20'
                    j = find(strcmp('wpo_n20',names));
                case 'wpon25'
                    % Western Pacific cruise not completely included,
                    % choose neighboring station instead.
                    j = find(strcmp('wpo_n30',names));
                case 'wpon30'
                    j = find(strcmp('wpo_n30',names));
                case 'wpos05'
                    % Western Pacific cruise not completely included,
                    % choose neighboring station instead.
                    j = find(strcmp('wpo_s10',names));
                case 'wpos10'
                    j = find(strcmp('wpo_s10',names));
                case 'wpos15'
                    % Western Pacific cruise not completely included,
                    % choose neighboring station instead.
                    j = find(strcmp('wpo_s20',names));
                case 'wpos20'
                    j = find(strcmp('wpo_s20',names));            
            end
            if ~isempty(j)
                LSCErespindex(i) = j;
            else
                error(['Could not find station: ' dname])
            end
            
        end
    end
    
    tsmon = tsmon(:,:,:,LSCErespindex);
    tspresub = tspresub(:,:,LSCErespindex);
    
end


% if the responses have 350 in (ie. offset ppmv) - subtract it off
if tsmon(1,1,1,1) > 300
    tsmon    = tsmon    - 350.0;
    tspresub = tspresub - 350.0;
end
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%