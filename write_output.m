function outcfg = write_output(invcfg,obscfg,outcfg,sources,data,toggle)  %#ok<INUSL>
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Writes the output for the NZ inversion.
%
% Author: Kay Steinkamp
% Date: Jan 2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define and check some parameters
% ----------------------------------------------------------------------
global monthnm weeknm

header = [invcfg.type ' NZ inversion with daily data'];
weekly = strcmp(invcfg.type,'weekly');
pstr = [num2str(invcfg.period(1)) '-' num2str(invcfg.period(2))];
srcnm = sources.name;
datnm([1 3 5]) = data.sitename;
datnm([2 4 6]) = {'sd'};
nsites = data.nsites;
nreg = sources.nreg;
noff = data.nsites;
srcunit = sources.unit;
day = cellstr(datestr(data.day));
nday = length(day);
ntime = data.ntimes;
for i = 1:ntime
    time{i} = data.type((18:25)+(i-1)*10);
end
datunit = data.unitsX;
if weekly
    nspatial = sources.nreg * sources.nweek;
else
    nspatial = sources.nreg * sources.nmon;
end
nt = sources.nt;
tname = pstr;
for y = 1:length(invcfg.period(1):invcfg.period(2))
    if weekly
        tname = [tname weeknm];
    else
        tname = [tname monthnm];
    end
end

if length(data.sitename) ~= nsites
    error('write_output: incompatible site/constraint numbers!')
end
if nspatial ~= sources.nsrc
    error('write_output: incompatible source numbers!')
end


% Output directory
% ----------------------------------------------------------------------
switch toggle  
    case 'prior'
        if ~exist([invcfg.workdir outcfg.writedir],'dir')
            mkdir([invcfg.workdir outcfg.writedir])
        end
        % Create run-specific output directory
        % (make sure no previous results get overwritten)
        runD = strcat([invcfg.workdir outcfg.writedir],date);
        runoutdir = [runD '/'];
        i = 1;
        while exist(runoutdir,'dir')
            runoutdir = [runD '.' int2str(i) '/'];
            i = i + 1;
        end
        mkdir(runoutdir)  % run-specific output folder
        % add to outcfg structure
        outcfg.runwritedir = runoutdir;
    case 'post'
        % Get output directory (known from prior function call)
        runoutdir = outcfg.runwritedir;
    otherwise
        error('write_output: toggle must read either "prior" or "post"!')
end


% Write source output
% ----------------------------------------------------------------------
if outcfg.dosources 
    
    switch toggle
        case 'prior'
            src = sources.prior;
            srcunc = sources.priorunc;
            srcuncPmean = sqrt(sum(srcunc.^2))/nt; % no correlations
            off = sources.priorV(nspatial+1:end);
            offunc = sources.prioruncV(nspatial+1:end);
            
        case 'post'
            src = sources.post;
            srcunc = sources.postunc;
            for r = 1:nreg
                rI = (r-1)*nt+1:r*nt;
                srcuncPmean(r) = mean(mean(sources.postcov(rI,rI)))^0.5;
            end
            off = sources.postV(nspatial+1:end);
            offunc = sources.postuncV(nspatial+1:end);
    end
    
    filename = ['outsrc_',toggle,'_',pstr,'.txt'];
    f15 = fopen([runoutdir filename],'wt');
    
    title1 = [toggle,' sources for period ',pstr];
    titleunit = ['regional sources in ' srcunit ', offset in ppm'];
    title2 = [toggle,' source uncertainty for period ',pstr];
    
    form1 = strcat('%19s%11s',strgrow('%8s',nt),'\n');
    form2 = strcat('%19s%11.2f',strgrow('%8.2f',nt),'\n');
    
    % write sources
    fprintf(f15,'%s\n',header);
    fprintf(f15,'\n');
    fprintf(f15,'%s\n',title1);
    fprintf(f15,'%s\n',titleunit);
    fprintf(f15,'%i,  %i\n',[nreg nt]);
    fprintf(f15,form1,'',tname{:});
    for i = 1:nreg
        fprintf(f15,form2,srcnm{i},mean(src(:,i)),src(:,i)');
    end
    % offset
    for i = 1:noff
        fprintf(f15,'%19s%11.2f\n',srcnm{nreg+i},off(i));
    end
    fprintf(f15,'\n');
    
    % write source uncertainty
    fprintf(f15,'%s\n',title2);
    fprintf(f15,'%s\n',titleunit);
    fprintf(f15,'%i,  %i\n',[nreg nt]);
    fprintf(f15,form1,'',tname{:});
    for i = 1:nreg
        fprintf(f15,form2,srcnm{i}, srcuncPmean(i), srcunc(:,i)');
    end
    % offset
    for i = 1:noff
        fprintf(f15,'%19s%11.2f\n',srcnm{nreg+i},offunc(i));
    end
    
    fclose(f15);
            
end


% Write data output
% ----------------------------------------------------------------------
if outcfg.dodata  
    
    switch toggle
        case 'prior'
            dat = data.priorX;
            datunc = data.priorXunc;
            
        case 'post'
            dat = data.postX;
            datunc = data.postXunc;
    end  
    
    filename = strcat('outdat_',toggle,'_',pstr,'.txt');
    f25 = fopen([runoutdir filename],'wt');
    
    title1 = [toggle,' data (deviations from baseline) and '...
        'standard deviations for period',' ',pstr];
    titleunit = ['station data are in ' datunit];
    
    form1 = strcat('%11s','%12s',strgrow('%10s',2*nsites),'\n');
    form2 = strcat('%11s','%12s',strgrow('%10.2f',2*nsites),'\n');
    
    % write data with standard deviations
    fprintf(f25,'%s\n',header);
    fprintf(f25,'\n');
    fprintf(f25,'%s\n',title1);
    fprintf(f25,'%s\n',titleunit);
    fprintf(f25,'%i,  %i\n',[nday*ntime 2*nsites]);
    fprintf(f25,form1,'Day','Starthour',datnm{:});
    datblock = [];
    for s = 1:nsites
        datblock = [datblock dat(:,s) datunc(:,s)];
    end
    for i = 1:nday
        for t = 1:ntime
            fprintf(f25,form2,day{i},time{t},datblock(t+(i-1)*ntime,:));
        end
    end
    fprintf(f25,'\n');
    
    fclose(f25);
end


% Save PDFs in mat files for use in further post-processing scripts 
% (e.g. make groups, mapping of fluxes, etc.)
% ----------------------------------------------------------------------
if strcmp(toggle,'post')
    
    filename = [runoutdir,'cfg_src_dat_',pstr,'.mat'];

    save(filename,'invcfg','obscfg','outcfg','sources','data')

end


% Move files written by earlier routines (called by
% 'regNZinv', eg. 'read_config') to appropriate directories
% ----------------------------------------------------------------------
if strcmp(toggle,'prior')   % needs to be done only once
    
    % configuration file
    file = strcat([invcfg.workdir outcfg.writedir],'configuration');
    file_new = strcat(runoutdir,'configuration');
    if exist(file,'file') && ~exist(file_new,'file')
        movefile(file, file_new)
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%