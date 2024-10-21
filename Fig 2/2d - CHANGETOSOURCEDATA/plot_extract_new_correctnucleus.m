clear all

global umperpix spf size_x size_y

umperpix = 354.25/512;
size_x = 512;size_y = 512;

%%

% figure;

for lighton = -1:1:1

    switch lighton

        case 0
            t_total = 50;
            spf = 2;
            taxis = (-t_total:-1).*spf;
            filename = ['slide2_pos3_tseries1_40x_mCherry_t50_spf_2s_lightoff.czi - C=0'];
            loaddir = ['save_kymocurvextract_' filename '_20221001 15_13'];

            load([loaddir '/' filename '_mem_cyto_pxlist.mat'],'mem_pxlist','cyto_pxlist',...
                'mem_pxlist_mean','mem_pxlist_std','cyto_pxlist_mean','cyto_pxlist_std');

            mem_pxlist_gef_t1 = mem_pxlist;
            cyto_pxlist_gef_t1 = cyto_pxlist;

            mem_pxlist_mean_gef_t1 = mem_pxlist_mean;
            cyto_pxlist_mean_gef_t1 = cyto_pxlist_mean;

            mem_pxlist_std_gef_t1 = mem_pxlist_std;
            cyto_pxlist_std_gef_t1 = cyto_pxlist_std;

            mem_pxlist_se_gef_t1 = mem_pxlist_std./sqrt(cellfun(@length,mem_pxlist));
            cyto_pxlist_se_gef_t1 = cyto_pxlist_std./sqrt(cellfun(@length,cyto_pxlist));

        case 1
            t_total = 376;
            spf = 9;
            taxis = (0:(t_total-1)).*spf;
            filename = ['slide2_pos3_tseries2_40x_mCherry_gfp_t376_spf_9s_lighton_viaimaging.czi - C=0'];
            loaddir = ['save_kymocurvextract_' filename '_20221001 15_28'];

            load([loaddir '/' filename '_mem_cyto_pxlist.mat'],'mem_pxlist','cyto_pxlist',...
                'mem_pxlist_mean','mem_pxlist_std','cyto_pxlist_mean','cyto_pxlist_std');

            mem_pxlist_gef_t2 = mem_pxlist;
            cyto_pxlist_gef_t2 = cyto_pxlist;

            mem_pxlist_mean_gef_t2 = mem_pxlist_mean;
            cyto_pxlist_mean_gef_t2 = cyto_pxlist_mean;

            mem_pxlist_std_gef_t2 = mem_pxlist_std;
            cyto_pxlist_std_gef_t2 = cyto_pxlist_std;

            mem_pxlist_se_gef_t2 = mem_pxlist_std./sqrt(cellfun(@length,mem_pxlist));
            cyto_pxlist_se_gef_t2 = cyto_pxlist_std./sqrt(cellfun(@length,cyto_pxlist));

        case -1 % Rho channel
            t_total = 376;
            spf = 9;
            taxis = (0:(t_total-1)).*spf;
            filename = ['slide2_pos3_tseries2_40x_mCherry_gfp_t376_spf_9s_lighton_viaimaging.czi - C=2'];
            loaddir = ['save_kymocurvextract_' filename '_20221001 16_13'];

            load([loaddir '/' filename '_mem_cyto_pxlist.mat'],'mem_pxlist','cyto_pxlist',...
                'mem_pxlist_mean','mem_pxlist_std','cyto_pxlist_mean','cyto_pxlist_std');

            mem_pxlist_rho_t2 = mem_pxlist;
            cyto_pxlist_rho_t2 = cyto_pxlist;

            mem_pxlist_mean_rho_t2 = mem_pxlist_mean;
            cyto_pxlist_mean_rho_t2 = cyto_pxlist_mean;

            mem_pxlist_std_rho_t2 = mem_pxlist_std;
            cyto_pxlist_std_rho_t2 = cyto_pxlist_std;

            mem_pxlist_se_rho_t2 = mem_pxlist_std./sqrt(cellfun(@length,mem_pxlist));
            cyto_pxlist_se_rho_t2 = cyto_pxlist_std./sqrt(cellfun(@length,cyto_pxlist));

    end
end

%%

figure;

t_total = 50;
spf = 2;
taxis_t1 = (-t_total:-1).*spf;
t_total = 376;
spf = 9;
taxis_t2 = (0:(t_total-1)).*spf;
taxis_full = [taxis_t1 taxis_t2];

fill([taxis_full flip(taxis_full)],...
    [[mem_pxlist_mean_gef_t1,mem_pxlist_mean_gef_t2]-[mem_pxlist_se_gef_t1,mem_pxlist_se_gef_t2],...
    flip([mem_pxlist_mean_gef_t1,mem_pxlist_mean_gef_t2]+[mem_pxlist_se_gef_t1,mem_pxlist_se_gef_t2])],...
    'c','facealpha',0.4,'EdgeColor','none');
hold on;
h1 = plot(taxis_full,[mem_pxlist_mean_gef_t1,mem_pxlist_mean_gef_t2]);hold on;

fill([taxis_full flip(taxis_full)],...
    [[cyto_pxlist_mean_gef_t1,cyto_pxlist_mean_gef_t2]-[cyto_pxlist_se_gef_t1,cyto_pxlist_se_gef_t2],...
    flip([cyto_pxlist_mean_gef_t1,cyto_pxlist_mean_gef_t2]+[cyto_pxlist_se_gef_t1,cyto_pxlist_se_gef_t2])],...
    'r','facealpha',0.4,'EdgeColor','none');
hold on;
h2 = plot(taxis_full,[cyto_pxlist_mean_gef_t1,cyto_pxlist_mean_gef_t2]);hold on;

legend([h1,h2],{'mem_gef','cyto_gef'},'location','southeast','fontsize',14);

xlim([-100 4400]);

% fig_current = gcf; fig_current.Renderer = 'painters';
% print(fig_current,['S4b_gef_dynamics_unnormalized'],'-dpdf');

%%

figure;

fill([taxis_t2 flip(taxis_t2)],...
    [mem_pxlist_mean_rho_t2-mem_pxlist_se_rho_t2,flip(mem_pxlist_mean_rho_t2+mem_pxlist_se_rho_t2)],...
    'c','facealpha',0.4,'EdgeColor','none');
hold on;
h1 = plot(taxis_t2,mem_pxlist_mean_rho_t2);hold on;

fill([taxis_t2 flip(taxis_t2)],...
    [cyto_pxlist_mean_rho_t2-cyto_pxlist_se_rho_t2,flip(cyto_pxlist_mean_rho_t2+cyto_pxlist_se_rho_t2)],...
    'r','facealpha',0.4,'EdgeColor','none');
hold on;
h2 = plot(taxis_t2,cyto_pxlist_mean_rho_t2);hold on;

plot([0,0],[9,14],'k--');hold on;

legend([h1,h2],{'mem_rho','cyto_rho'},'location','southeast','fontsize',14);

xlim([-100 4400]);

% fig_current = gcf; fig_current.Renderer = 'painters';
% print(fig_current,['S4b_rho_dynamics_unnormalized'],'-dpdf');

%%

figure;
mem_tot_gef_t1 = zeros(1,length(taxis_t1));
cyto_tot_gef_t1 = zeros(1,length(taxis_t1));
for t = 1:1:length(taxis_t1)
    mem_tot_gef_t1(t) = sum(mem_pxlist_gef_t1{t});
    cyto_tot_gef_t1(t) = sum(cyto_pxlist_gef_t1{t});
end
mem_tot_gef_t2 = zeros(1,length(taxis_t2));
cyto_tot_gef_t2 = zeros(1,length(taxis_t2));
mem_tot_rho_t2 = zeros(1,length(taxis_t2));
cyto_tot_rho_t2 = zeros(1,length(taxis_t2));
for t = 1:1:length(taxis_t2)
    mem_tot_gef_t2(t) = sum(mem_pxlist_gef_t2{t});
    cyto_tot_gef_t2(t) = sum(cyto_pxlist_gef_t2{t});
    mem_tot_rho_t2(t) = sum(mem_pxlist_rho_t2{t});
    cyto_tot_rho_t2(t) = sum(cyto_pxlist_rho_t2{t});
end

array_gef = ([mem_tot_gef_t1,mem_tot_gef_t2])./...
    (([mem_tot_gef_t1,mem_tot_gef_t2])+([cyto_tot_gef_t1,cyto_tot_gef_t2]));
h1 = plot(taxis_full,(array_gef-mean(array_gef(1:50)))./(max(array_gef)-mean(array_gef(1:50))));hold on;
array_rho = mem_tot_rho_t2./(mem_tot_rho_t2+cyto_tot_rho_t2);
h2 = plot(taxis_t2,(array_rho-array_rho(1))./(max(array_rho)-array_rho(1)));hold on;

legend([h1,h2],{'mem_frac_gef','mem_frac_rho'},'location','southeast','fontsize',14);

xlim([-100 4400]);ylim([0 1]);

fig_current = gcf; fig_current.Renderer = 'painters';
print(fig_current,['1gnew_gef_rho_dynamics_normalized_20230215'],'-dpdf');


