
filenamelist_wt = {'d7oo11wt','d7oo3wt','d7oo5wt','d7oo8wt','d7oo7wt'};
speeds_wt = zeros(1,length(filenamelist_wt));
for k = 1:1:5
    filename = filenamelist_wt{k};
    load(['wavespeed final wt/' filename '.mat'],'speed');
    speeds_wt(k) = speed;
end

filenamelist_larg = {'dxoo4', 'd4oo1_t1', 'd4oo1_t2', 'd0oo2_t1', 'd0oo2_t2'};
speeds_larg = zeros(1,length(filenamelist_larg));
for k = 1:1:5
    filename = filenamelist_larg{k};
    load(['wavespeed final larg/' filename '.mat'],'speed');
    speeds_larg(k) = speed;
end

filenamelist_ect2 = {'d0oo1','d1oo1','d2oo1','d3oo1'};
speeds_ect2 = zeros(1,length(filenamelist_ect2));
for k = 1:1:4
    filename = filenamelist_ect2{k};
    load(['wavespeed final ect2/' filename '.mat'],'speed');
    speeds_ect2(k) = speed;
end

slows = [0.079781;0.0796885;0.0796111;0.0795419;0.0794734;0.0794136;0.079355;0.0792894;0.0792331;0.0791782;0.0791248;0.079073;0.079023;0.0789747;0.0789285;0.0788843;0.0788423;0.0788027;0.0787559;0.0787114;0.0786693;0.0786299;0.0785932;0.0785596;0.0785184;0.07848;0.0784449;0.0784131;0.0783848;0.0783487;0.078316;0.0782872;0.0782501;0.0782167;0.0781873;0.0781492;0.0781151;0.0780853;0.0780602;0.0780401;0.0780109;0.0779867;0.0779528;0.077924;0.0779006;0.0778671;0.077839;0.0778168;0.077801;0.0777921;0.0777549;0.0777238;0.0776992;0.0776817;0.0776719;0.0776505;0.0776165;0.0775899;0.0775712;0.0775612;0.0775383;0.0775242;0.0774962;0.0774773;0.0774684;0.0774448;0.0774315;0.0774026;0.0773841;0.0773771;0.0773537;0.0773422;0.0773438;0.0773282;0.0772935;0.0772715;0.0772636;0.0772713;0.0772595;0.0772259;0.0772076;0.0772061;0.0772238;0.0772193;0.0771899;0.0771791;0.0771893;0.0771733;0.0771276;0.077102;0.0770996;0.0771237;0.0771173;0.077076;0.0770605;0.0770747;0.0771234;0.0770622;0.0770297;0.077031;0.077072;0.0770693;0.0770148;0.076999;0.0770294;0.0771151;0.0770323;0.0769941;0.0770107;0.0770942;0.0771108;0.0770441;0.0770436;0.0771263;0.0771189;0.0769964;0.076953;0.0770132;0.0772091;0.0770028;0.0768922;0.0769144;0.0771201;0.0771529;0.0769369;0.0768984;0.0771261;0.0770634;0.0765463;0.0762393;0.0763043;0.0770186;0.0772932;0.076677;0.0767647;0.0784615;0.0793706];
wts = [0.52005;0.499883;0.48179;0.471048;0.463405;0.457742;0.452677;0.448633;0.445033;0.442247;0.439272;0.436137;0.434;0.431169;0.428872;0.427864;0.427206;0.427473;0.425524;0.424242;0.42619;0.422857;0.43];
fasts = [3.75167;3.55357;3.08];

slows = slows.*1.2*60; % convert to micron per min and convert diameter 100 to 120
wts = wts.*1.2*60;
fasts = fasts.*1.2*60;

figure;
scatter(1,mean(speeds_wt),100,'ks');hold on;
errorbar(1,mean(speeds_wt),std(speeds_wt,1),std(speeds_wt,1),'k-');hold on;
scatter(1,mean(wts),100,'ks','filled');hold on;
errorbar(1,mean(wts),std(wts,1),std(wts,1),'k-');hold on;

scatter(2,mean(speeds_larg),100,'ko');hold on;
errorbar(2,mean(speeds_larg),std(speeds_larg,1),std(speeds_larg,1),'k-');hold on;
scatter(2,mean(slows),100,'ko','filled');hold on;
errorbar(2,mean(slows),std(slows,1),std(slows,1),'k-');hold on;

scatter(3,mean(speeds_ect2),100,'k^');hold on;
errorbar(3,mean(speeds_ect2),std(speeds_ect2,1),std(speeds_ect2,1),'k-');hold on;
scatter(3,mean(fasts),100,'k^','filled');hold on;
errorbar(3,mean(fasts),std(fasts,1),std(fasts,1),'k-');hold on;

xlim([0.5 3.5]);ylim([4 400]);

box on

set(gca,'yscale','log');

fig_current = gcf; fig_current.Renderer = 'painters';
print(fig_current,['20230219_wave_speed_log'],'-dpdf');