function makeplots(invcfg, outcfg, data, sources)   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Selection of plots, to be called after inversion.
%
% Author: Kay Steinkamp
% Date: Feb 2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% close all

% check if statistics toolbox is available
Slic = license('checkout','statistics_toolbox');

%% plot regional annual sources
% ----------------------------------------------------------------------
figure('Position',[20 60 900 800],'color','w')
set(gcf,'renderer','zbuffer','PaperPositionMode','auto')
for y = 1:length(invcfg.period(1):invcfg.period(2))
    subplot(length(invcfg.period(1):invcfg.period(2)),1,y);
    
    h = bar([sources.priorA{y};sources.postA{y}]');
    title([num2str(invcfg.period(1)+y-1) ' annual regional sources'])
    UlinesPrior = [sources.priorA{y} - sources.prioruncA{y};...
        sources.priorA{y} + sources.prioruncA{y}];
    UlinesPost = [sources.postA{y} - sources.postuncA{y};...
        sources.postA{y} + sources.postuncA{y}];

    % add the uncertainty lines
    for r = 1:sources.nreg;
        line([r r]-1/7,UlinesPrior(:,r),'LineStyle','-','Color',[0.5 0.5 0.5])
        line([r r]+1/7,UlinesPost(:,r),'LineStyle','-','Color',[0.5 0.5 0.5])
    end
    set(gca,'xlim',[0.5 sources.nreg+0.5])
    set(gca,'xtick',0.5:sources.nreg+0.5,'xticklabel',[])
    grid on; axis tight
    yl = get(gca,'ylim');
    ypos = yl(1)-0.05*diff(yl);
    text((1:sources.nreg)-1/4,ypos*ones(sources.nreg,1),num2cell(1:sources.nreg))
    ylabel('Mt CO2 yr-1')
    legend(h,'prior','posterior')
end
try
    print('-dpng',[outcfg.runwritedir 'regional_annual_sources.png'])
catch %#ok<*CTCH>
    disp('Could not save plot')
end


%% plot source update
% ----------------------------------------------------------------------
if Slic
    figure('Position',[20 60 1000 600],'color','w')
    set(gcf,'renderer','zbuffer','PaperPositionMode','auto')
    subplot(1,2,1)
    boxplot((sources.post-sources.prior),'plotstyle','compact','symbol','b')
    % set(gca,'ylim',[0 101])
    % ylabel('% of prior source')
    ylabel('Mt CO2 yr-1')
    xlabel('Region No.')
    title('Source update, all sources in 2011-2012')
    yl = get(gca,'ylim');
    bins = linspace(yl(1),yl(2),30);
    subplot(1,2,2)
    nn = histc(sources.post(:)-sources.prior(:),bins);
    barh(bins,nn)
    try
        print('-dpng',[outcfg.runwritedir 'regional_source_update.png'])
    catch
        disp('Could not save plot')
        
    end
end


%% plot regional weekly sources 
% (if inversion type is 'weekly')
% ----------------------------------------------------------------------
if strcmp(invcfg.type,'weekly')
    pri = sources.prior;
    priU = sources.priorunc;
    pos = sources.post;
    posU = sources.postunc;
    figure('Position',[20 130 1450 850],'color','w')
    set(gcf,'renderer','zbuffer','PaperPositionMode','auto')
    for r = 1:sources.nreg
        rfac = 1;
        unitstr = 'Mt CO2 yr-1';
%         rfac = 1e9/sources.regarea(r);
%         unitstr = 'kg CO2 m-2 yr-1';
        subplot(4,7,r);
        if all(priU(:,r)<1000)
            patch([1:sources.nweek sources.nweek:-1:1]', ...
                rfac*[pri(:,r)-priU(:,r);...
                pri(end:-1:1,r)+priU(end:-1:1,r)],...
                'w','EdgeColor','none','FaceColor',[.8 .8 .8])
        end
        hold on
        patch([1:sources.nweek sources.nweek:-1:1]', ...
            rfac*[pos(:,r)-posU(:,r);...
            pos(end:-1:1,r)+posU(end:-1:1,r)],...
            'w','EdgeColor','none','FaceColor',[1 .7 .7])
        plot([1 sources.nweek],[0 0],'k:')
%         if r<16
%             set(gca,'ylim',[-10 10]);
%         else
%             set(gca,'ylim',[-3 3]);
%         end
           
        yl = get(gca,'ylim');
        for i = 1:length(sources.years)-1
            plot([1 1]*i*48.5,yl,'w')
        end
        plot(rfac*pri(:,r),'k')
        plot(rfac*pos(:,r),'r')
        plot(rfac*(pos(:,r)-pri(:,r)),'g')
       %         plot(smooth(1:sources.nweek,rfac*pos(:,r),0.3,'loess'),...
        %             'r','linewidth',2);
%         plot(smooth(rfac*pos(:,r),4),'r','linewidth',2);
        set(gca,'xlim',[1 sources.nweek],'xtick',[])
        xlabel('Week in period')
        ylabel(unitstr)
        title(sources.name{r})
    end
    
    % add NZ totals
    subplot(4,7,26);
    plot([1 sources.nweek],[0 0],'k:')
    hold on
    plot(sum(pri(:,1:8),2),'k')
    plot(sum(pos(:,1:8),2),'r')
%     plot(smooth(1:sources.nweek,sum(pos(:,1:8),2),0.3,'loess'),'r','linewidth',2);
    yl = get(gca,'ylim');
    plot([1 1]*(sources.nweek+1)/2,yl,'color',[.8 .8 .8])
    set(gca,'xlim',[1 sources.nweek],'xtick',[])
    xlabel('Week in period')
    ylabel('Mt CO2 yr-1')
    title('NZ North Island')
    text(2,yl(1)+0.09*diff(yl),[num2str(mean(sum(pri(:,1:8),2))) 'A-Prior'])
    text(2,yl(1)+0.02*diff(yl),[num2str(mean(sum(pos(:,1:8),2))) 'A-Post'])
    
    subplot(4,7,27);
    plot([1 sources.nweek],[0 0],'k:')
    hold on
    plot(sum(pri(:,9:15),2),'k')
    plot(sum(pos(:,9:15),2),'r')
%     plot(smooth(1:sources.nweek,sum(pos(:,9:15),2),0.3,'loess'),'r','linewidth',2);
    yl = get(gca,'ylim');
    plot([1 1]*(sources.nweek+1)/2,yl,'color',[.8 .8 .8])
    set(gca,'xlim',[1 sources.nweek],'xtick',[])
    xlabel('Week in period')
    ylabel('Mt CO2 yr-1')
    title('NZ South Island')
    text(2,yl(1)+0.09*diff(yl),[num2str(mean(sum(pri(:,9:15),2))) 'A-Prior'])
    text(2,yl(1)+0.02*diff(yl),[num2str(mean(sum(pos(:,9:15),2))) 'A-Post'])
    
    subplot(4,7,28);
    plot([1 sources.nweek],[0 0],'k:')
    hold on
    plot(sum(pri(:,1:15),2),'k')
    plot(sum(pos(:,1:15),2),'r')
%     plot(smooth(1:sources.nweek,sum(pos(:,1:15),2),0.3,'loess'),'r','linewidth',2);
    yl = get(gca,'ylim');
    plot([1 1]*(sources.nweek+1)/2,yl,'color',[.8 .8 .8])
    set(gca,'xlim',[1 sources.nweek],'xtick',[])
    xlabel('Week in period')
    ylabel('Mt CO2 yr-1')
    title('NZ Total')
    text(2,yl(1)+0.09*diff(yl),[num2str(mean(sum(pri(:,1:15),2))) 'A-Prior'])
    text(2,yl(1)+0.02*diff(yl),[num2str(mean(sum(pos(:,1:15),2))) 'A-Post'])
    try
        print('-dpng',[outcfg.runwritedir 'regional_weekly_sources.png'])
    catch
        disp('Could not save plot')
    end
end


%% plot regional monthly sources
% ----------------------------------------------------------------------
if strcmp(invcfg.type,'weekly')
    pri = sources.priorM;
    priU = sources.prioruncM;
    pos = sources.postM;
    posU = sources.postuncM;
else
    pri = sources.prior;
    priU = sources.priorunc;
    pos = sources.post;
    posU = sources.postunc;
end
figure('Position',[20 130 1450 850],'color','w')
set(gcf,'renderer','zbuffer','PaperPositionMode','auto')
for r = 1:sources.nreg
    rfac = 1;
    unitstr = 'Mt CO2 yr-1';
%     rfac = 1e9/sources.regarea(r);
%     unitstr = 'kg CO2 m-2 yr-1';
    subplot(4,7,r); 
    patch([1:sources.nmon sources.nmon:-1:1]', ...
        rfac*[pri(:,r)-priU(:,r);...
        pri(end:-1:1,r)+priU(end:-1:1,r)],...
        'w','EdgeColor','none','FaceColor',[.8 .8 .8])
    hold on
    patch([1:sources.nmon sources.nmon:-1:1]', ...
        rfac*[pos(:,r)-posU(:,r);...
        pos(end:-1:1,r)+posU(end:-1:1,r)],...
        'w','EdgeColor','none','FaceColor',[1 .7 .7])
    plot([1 sources.nmon],[0 0],'k:')
%     if r<16
%         set(gca,'ylim',[-8 8]);
%     else
%         set(gca,'ylim',[-3 3]);
%     end
    yl = get(gca,'ylim');
    for i = 1:length(sources.years)-1
        plot([1 1]*i*12.5,yl,'w')
    end
    plot(rfac*pri(:,r),'k')
    plot(rfac*pos(:,r),'r')
    plot(smooth(1:sources.nmon,rfac*pos(:,r),0.3,'loess'),...
        'r','linewidth',2);
%     plot(smooth(pos(:,r),3),'r','linewidth',2);
    set(gca,'xlim',[1 sources.nmon],'xtick',[])
    xlabel('Month in period')
    ylabel(unitstr)
    title(sources.name{r})
end

% add NZ totals
subplot(4,7,26);
plot([1 sources.nmon],[0 0],'k:')
hold on
plot(sum(pri(:,1:8),2),'k')
plot(sum(pos(:,1:8),2),'r')
plot(smooth(1:sources.nmon,sum(pos(:,1:8),2),0.3,'loess'),'r','linewidth',2);
yl = get(gca,'ylim');
plot([1 1]*(sources.nmon+1)/2,yl,'color',[.8 .8 .8])
set(gca,'xlim',[1 sources.nmon],'xtick',[])
xlabel('Month in period')
ylabel('Mt CO2 yr-1')
title('NZ North Island')
text(2,yl(1)+0.09*diff(yl),[num2str(mean(sum(pri(:,1:8),2))) 'A-Prior'])
text(2,yl(1)+0.02*diff(yl),[num2str(mean(sum(pos(:,1:8),2))) 'A-Post'])

subplot(4,7,27);
plot([1 sources.nmon],[0 0],'k:')
hold on
plot(sum(pri(:,9:15),2),'k')
plot(sum(pos(:,9:15),2),'r')
plot(smooth(1:sources.nmon,sum(pos(:,9:15),2),0.3,'loess'),'r','linewidth',2);
yl = get(gca,'ylim');
plot([1 1]*(sources.nmon+1)/2,yl,'color',[.8 .8 .8])
set(gca,'xlim',[1 sources.nmon],'xtick',[])
xlabel('Month in period')
ylabel('Mt CO2 yr-1')
title('NZ South Island')
text(2,yl(1)+0.09*diff(yl),[num2str(mean(sum(pri(:,9:15),2))) 'A-Prior'])
text(2,yl(1)+0.02*diff(yl),[num2str(mean(sum(pos(:,9:15),2))) 'A-Post'])

subplot(4,7,28);
plot([1 sources.nmon],[0 0],'k:')
hold on
plot(sum(pri(:,1:15),2),'k')
plot(sum(pos(:,1:15),2),'r')
plot(smooth(1:sources.nmon,sum(pos(:,1:15),2),0.3,'loess'),'r','linewidth',2);
yl = get(gca,'ylim');
plot([1 1]*(sources.nmon+1)/2,yl,'color',[.8 .8 .8])
set(gca,'xlim',[1 sources.nmon],'xtick',[])
xlabel('Month in period')
ylabel('Mt CO2 yr-1')
title('NZ Total')
text(2,yl(1)+0.09*diff(yl),[num2str(mean(sum(pri(:,1:15),2))) 'A-Prior'])
text(2,yl(1)+0.02*diff(yl),[num2str(mean(sum(pos(:,1:15),2))) 'A-Post'])
try
    print('-dpng',[outcfg.runwritedir 'regional_monthly_sources.png'])
catch
    disp('Could not save plot')
end


%% plot uncertainty reduction
% ----------------------------------------------------------------------
if Slic
    figure('Position',[20 60 900 800],'color','w')
    set(gcf,'renderer','zbuffer','PaperPositionMode','auto')
    boxplot((sources.priorunc./sources.priorunc)*100)
    hold on
    if strcmp(invcfg.type,'weekly')
        boxplot((sources.postuncM./sources.prioruncM)*100,'symbol','r',...
            'plotstyle','compact','colors','g') % monthly
    end
    set(gca,'ylim',[0 101],'xticklabel','none')
    boxplot((cell2mat(sources.postuncA')./cell2mat(sources.prioruncA'))*100,...
        'symbol','k','plotstyle','compact','colors','m') % annual
    set(gca,'xticklabel','none')
    boxplot((sources.postunc./sources.priorunc)*100,'symbol','r') % weekly(monthly)
    set(gca,'ylim',[0 101])
    % set(gca,'ylim',[0 101],'xticklabel',cellstr(num2str((1:sources.nreg)')))
    ylabel('% of prior uncertainty')
    xlabel('Region No.')
    title('Uncertainty reduction, all sources in period')
    try
        print('-dpng',[outcfg.runwritedir 'regional_unc_reduction.png'])
    catch
        disp('Could not save plot')
    end
end


%% plot the posterior source uncertainty distribution
% ----------------------------------------------------------------------
figure('Position',[100 650 550 500])
hist(sources.postuncV(1:sources.nsrc),30)
xlabel('Mt CO2 yr-1')
txtstr = ['median: ' num2str(median(sources.postuncV(1:sources.nsrc)))...
    ' Mt CO2 yr-1'];
text(mean(get(gca,'xtick')),mean(get(gca,'ytick')),txtstr)
txtstr = ['mean: ' num2str(mean(sources.postuncV(1:sources.nsrc)))...
    ' Mt CO2 yr-1'];
text(mean(get(gca,'xtick')),mean(get(gca,'ytick'))-10,txtstr)
try
    print('-dpng',[outcfg.runwritedir 'histogram_posterior_Sunc.png'])
catch
    disp('Could not save plot')
end

if strcmp(invcfg.type,'weekly')
    figure('Position',[120 650 550 500])
    hist(sources.postuncM(:),30)
    xlabel('Mt CO2 yr-1')
    txtstr = ['median: ' num2str(median(sources.postuncM(:))) ' Mt CO2 yr-1'];
    txtstr2 = ['median w/o error corr: ' ...
        num2str(median(sources.postuncV(1:sources.nsrc)/sqrt(4))) ' Mt CO2 yr-1'];
    text(mean(get(gca,'xtick'))-15,mean(get(gca,'ytick')),txtstr)
    text(mean(get(gca,'xtick'))-15,mean(get(gca,'ytick'))-10,txtstr2)
    txtstr = ['mean: ' num2str(mean(sources.postuncM(:))) ' Mt CO2 yr-1'];
    txtstr2 = ['mean w/o error corr: ' ...
        num2str(mean(sources.postuncV(1:sources.nsrc)/sqrt(4))) ' Mt CO2 yr-1'];
    text(mean(get(gca,'xtick'))-15,mean(get(gca,'ytick'))-20,txtstr)
    text(mean(get(gca,'xtick'))-15,mean(get(gca,'ytick'))-30,txtstr2)
    try
        print('-dpng',[outcfg.runwritedir 'histogram_posterior_Sunc_monthly.png'])
    catch
        disp('Could not save plot')
    end
end


%% plot the posterior error correlation matrix
% (only monthly - weekly matrix is big)
% ----------------------------------------------------------------------
if strcmp(invcfg.type,'weekly')
    cov = sources.postcovM;
    N = sources.nmon*sources.nreg;
else
    cov = sources.postcov;
    N = sources.nsrc;
end
figure('Position',[100 650 550 500],'renderer','zbuffer')
corrmat = cov2corr_own(cov(1:N,1:N)) - eye(N);
pcolor_own(1:N,1:N,corrmat)
shading flat
colormap(own_colormaps(4))
set(gca,'clim',[-1 1]*max(abs(corrmat(:))))
colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%