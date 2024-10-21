
figure;

load('1C.mat','raw_memfrac','spf','t_on');

for k = 1:1:2

    tnum_expt = length(raw_memfrac{k});
    taxis_expt = ((0:1:tnum_expt-1)-t_on).*spf;

    tseries_expt_norm = (raw_memfrac{k}-mean(raw_memfrac{k}(1:(t_on-1))))./(max([raw_memfrac{:}])-mean(raw_memfrac{k}(1:(t_on-1))));

    plot(taxis_expt,tseries_expt_norm,'linewidth',1.5);hold on;
end

filename = 'photosensitvie_tag_attachment_cleaned';

load([filename '.mat'],'datatable','taxis');

for k = 1:1:2

    tseries_simu_norm  = (datatable{k}-datatable{k}(1))./(max([datatable{:}])-datatable{k}(1));

    plot(taxis,tseries_simu_norm,'linewidth',1.5);hold on;
end

xlim([-100 4400]);ylim([0,1]);

set(gca,'xtick',[0,1000,2000,3000,4000],'ytick',[0,0.2,0.4,0.6,0.8,1]);

pbaspect([1.2 1 1]);

plot([0,0],[0,1],'k--');hold on;
plot([4000,4000],[0,1],'k--');hold on;

fig_current = gcf; fig_current.Renderer = 'painters';
print(fig_current,['2C_activation_kinetics_global_expt_simu'],'-dpdf');

