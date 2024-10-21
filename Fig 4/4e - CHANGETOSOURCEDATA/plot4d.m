
filenamelist = {'dxoo6','dxoo3','dxoo4','d5oo1_t1','d5oo1_t3','d5oo1_t5','d4oo1_t1','d4oo1_t2','d0oo2_t2','d0oo2_tlong'}; % 'd0oo2_t1',
pinchorwavelist = [0,0,1,0,0,0,1,1,1,1];

figure;
for k_file = 1:1:10

    filename = filenamelist{k_file};

    %     load(['S6g plot/' filename '.mat'],'taxis','phiminlist','spf','umperpix','perimeter','t_off_local',...
    %         'phiplot','phiplot_lower','phiplot_upper','idx','phidot','phidot_errlh','phidot_errrh');

    load(['S6g plot/' filename '_prodist0.mat'],'progdist0','err_progdist0');
    load(['/Users/jinghui/Documents/Lab/light activation/2023Jan_SIFigs/S6/S6f density recruitment/optofolds final/' filename '.mat'],'optofolds');

    %     subplot(1,2,1);
    %     scatter(optofolds,phidot,80,'k');hold on;
    %     errorbar(optofolds,phidot,(phidot-phidot_errrh),(phidot_errlh-phidot),0,0,'linewidth',1.5,'color','k');hold on;
    %     xlim([1,3]);ylim([0 2*pi]);box on;
    %     ylabel('Angle (radian)','fontsize',18);
    %     xlabel('Opto-GEF* fold increase','fontsize',18);
    %     pbaspect([1,1,1]);
    %
    %     subplot(1,2,2);
    if pinchorwavelist(k_file)
        mycolor = 'r';
    else
        mycolor = 'k';
    end
    scatter(mean(optofolds),progdist0,100,mycolor);hold on;
    errorbar(mean(optofolds),progdist0,err_progdist0,err_progdist0,...
        std(optofolds,1),std(optofolds,1),'linewidth',1.5,'color',mycolor);hold on;


end

xlim([1,3]);ylim([0 650]);
box on;
ylabel('Distance (um)','fontsize',18);
xlabel('Opto-GEF* fold increase','fontsize',18);
pbaspect([1,1,1]);

load('/Users/jinghui/Documents/Lab/light activation/2023Jan_SIFigs/S6/S6g f theory/data.mat','xlist_rescaled','ylist_rescaled');
plot(xlist_rescaled,ylist_rescaled,'linewidth',1.5);hold on;
% plot(xlist_rescaled(xlist_rescaled<2.5),ylist_rescaled(xlist_rescaled<2.5),'linewidth',1.5);hold on;

fig_current = gcf; fig_current.Renderer = 'painters';
print(fig_current,['4d'],'-dpdf');
