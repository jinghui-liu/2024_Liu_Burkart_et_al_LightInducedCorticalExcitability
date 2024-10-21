%% check bws from raw dic video

clear all

global spf umperpix Lmax size_col size_row

filename = 'optoect220220307';

umperpix = 266.74/512;

graylim = [5 40];

size_col = 512; size_row = 512;

Lmax = 150/umperpix;

t_total = 1131;

t_on = 101;
t_off = 572;

spf = 8;

taxis = (0:1:(t_total-1)).*spf/60;

% dictifpath = '/Volumes/One Touch/2022 Photoactivation - Analysis/Analysis_WildtypeSCW/d7oo1-14wt/d7oo11wt_dic8bit.tif';
geftifpath = '/Volumes/One Touch/2022 Photoactivation - Analysis/2 Analysis_ReversibleContractions/optoect220220307/d2oo2_optoect220220307.tif';

% load(['bws_updated/bws_' filename '.mat'],'bws_x','bws_y','rad','thresh','smfac','sFac');

load('save_kymocurvextract_optoect220220307_20220911 16_30/boundary.mat','bws');

bws_x = cell(1,t_total);bws_y = cell(1,t_total);
for t = 1:1:t_total
    bws_x{t} = bws{t}(:,2);
    bws_y{t} = bws{t}(:,1);
end

ROI1 = [36,96;222,261];
ROI2 = [436,496;222,262];

%% change to polar coordinate and visualize

load('save_kymocurvextract_optoect220220307_20220911 16_30/kymos.mat','kymos');

% figure;imagesc(kymo);colorbar;

[xCoM_raw,yCoM_raw] = getxyCoM(bws_x,bws_y); % in pixels

t = 1:1:t_total;

drawon = 0;
[rhos_raw,phis_raw,I_shift] = getpolarCoM(t,bws_x,bws_y,xCoM_raw,yCoM_raw,drawon); % convert to polar coordinate & in microns

kymos_shifted = nan(size(kymos));
for t = 1:1:t_total
    N_bd = length(bws_x{t});
    kymo = kymos(1:N_bd,t);
    
    kymos_shifted(1:N_bd,t) = kymo(I_shift{t});
end



%% compare interpolate using interp1 (interp twice on 180 degree apart breakings on the circle and combine)

N_phiinterp = 3600; % 0.1 degree
phis_interpcor_minuspitopi = linspace(-pi,pi,N_phiinterp+1);
phis_interpcor_0to2pi = linspace(0,2*pi,N_phiinterp+1);

kymos_interp_minuspitopi = zeros(N_phiinterp+1,t_total); 
kymos_interp_0to2pi = zeros(N_phiinterp+1,t_total);

rhos_interp_minuspitopi = zeros(N_phiinterp+1,t_total); 
rhos_interp_0to2pi = zeros(N_phiinterp+1,t_total);

for t = 1:1:t_total
    
    N_bd = length(bws_x{t});
    
    phit = phis_raw{t};
    [uniq_phit,I_phi] = unique(phit);
    kymo = kymos_shifted(1:N_bd,t);
    kymo_interp = interp1(uniq_phit,kymo(I_phi),phis_interpcor_minuspitopi);
    kymos_interp_minuspitopi(:,t) = kymo_interp;
    
    rhot = rhos_raw{t}; % in microns
    rhot_interp = interp1(uniq_phit,rhot(I_phi),phis_interpcor_minuspitopi);
    rhos_interp_minuspitopi(:,t) = rhot_interp;
    
    phit2 = phis_raw{t};
    phit2(phit2<0) = phit2(phit2<0)+2*pi;
    [uniq_phit,I_phi] = unique(phit2);
    kymo = kymos_shifted(1:N_bd,t);
    kymo_interp = interp1(uniq_phit,kymo(I_phi),phis_interpcor_0to2pi);
    kymos_interp_0to2pi(:,t) = kymo_interp;
    
    rhot2 = rhos_raw{t}; % in microns
    rhot_interp = interp1(uniq_phit,rhot(I_phi),phis_interpcor_0to2pi);
    rhos_interp_0to2pi(:,t) = rhot_interp;
    
end

phis_interpcor= linspace(-pi,pi,N_phiinterp+1);

kymos_interp = nan(N_phiinterp+1,t_total);% on -\pi to \pi
kymos_interp1 = kymos_interp_minuspitopi;
kymos_interp2 = kymos_interp_0to2pi;kymos_interp2 = [kymos_interp2(N_phiinterp/2+1:N_phiinterp+1,:);kymos_interp2(1:N_phiinterp/2,:)];

mask_kymo1nan = (isnan(kymos_interp1))&(~isnan(kymos_interp2));
mask_kymo2nan = (~isnan(kymos_interp1))&(isnan(kymos_interp2));
mask_nonan = (~isnan(kymos_interp1))&(~isnan(kymos_interp2));
kymos_interp(mask_kymo1nan) = kymos_interp2(mask_kymo1nan);
kymos_interp(mask_kymo2nan) = kymos_interp1(mask_kymo2nan);
kymos_interp(mask_nonan) = (kymos_interp1(mask_nonan)+kymos_interp1(mask_nonan)).*0.5;

rhos_interp = nan(N_phiinterp+1,t_total);% on -\pi to \pi
rhos_interp1 = rhos_interp_minuspitopi;
rhos_interp2 = rhos_interp_0to2pi;rhos_interp2 = [rhos_interp2(N_phiinterp/2+1:N_phiinterp+1,:);rhos_interp2(1:N_phiinterp/2,:)];

mask_rho1nan = (isnan(rhos_interp1))&(~isnan(rhos_interp2));
mask_rho2nan = (~isnan(rhos_interp1))&(isnan(rhos_interp2));
mask_nonan = (~isnan(rhos_interp1))&(~isnan(rhos_interp2));
rhos_interp(mask_rho1nan) = rhos_interp2(mask_rho1nan);
rhos_interp(mask_rho2nan) = rhos_interp1(mask_rho2nan);
rhos_interp(mask_nonan) = (rhos_interp1(mask_nonan)+rhos_interp1(mask_nonan)).*0.5;

% figure;
% imagesc(kymos_interp_minuspitopi);colorbar;pbaspect([2,1,1]);hold on;
% title('kymo interp -\pi to \pi','fontsize',16);
% figure;
% imagesc(kymos_interp_0to2pi);colorbar;pbaspect([2,1,1]);hold on;
% title('kymo interp 0 to 2\pi','fontsize',16);
% figure;
% imagesc(kymos_interp);colorbar;pbaspect([2,1,1]);hold on;
% title('kymo interp combined (-\pi to \pi)','fontsize',16);

%% VP1, kymo local 

VPangle = 108-180;% in degree
shiftVP = findequal(phis_interpcor,VPangle*pi/180,1);

kymoss_interp = [kymos_interp;kymos_interp;kymos_interp];
kymos_interp_shiftVP = kymoss_interp(N_phiinterp+shiftVP-N_phiinterp/2:N_phiinterp+shiftVP+N_phiinterp/2,:);

figure;
imagesc(taxis,phis_interpcor,kymos_interp_shiftVP);colorbar;pbaspect([2,1,1]);hold on;
% 
% plot(taxis,(0+max([VPangle1_upper,VPangle1_lower])-mean(VPangle1)).*ones(1,length(taxis)),'w--');hold on;
% plot(taxis,(0-min([VPangle1_upper,VPangle1_lower])-mean(VPangle1)).*ones(1,length(taxis)),'w--');hold on;
% 
set(gca,'ytick',[-pi,0,pi],'yticklabel',{'-\pi','0','\pi'},'fontsize',14);
title(['kymo ' filename],'fontsize',14);
xlabel('Time (min)','fontsize',16); 
ylabel('Polar angle \phi','fontsize',16);

%%

kymos_interp_shiftVP0 = mean(kymos_interp_shiftVP(:,1:60),2)*ones(1,t_total);

% kymos_interp_shiftVP0 = kymos_interp_shiftVP; %[-(kymos_interp_shiftVP(:,74)-kymos_interp_shiftVP(:,73))*ones(1,73) zeros(size(kymos_interp_shiftVP,1),size(kymos_interp_shiftVP,2)-73)];
% for i = 1:1:size(kymos_interp_shiftVP0,2)
%     kymos_interp_shiftVP(:,i) = kymos_interp_shiftVP0(:,i)./mean(kymos_interp_shiftVP0(:,i));
% end


figure;

imagesc(taxis,phis_interpcor,kymos_interp_shiftVP-kymos_interp_shiftVP0);colorbar;hold on;
% plot(taxis,(0+max([VPangle1_upper,VPangle1_lower])-mean(VPangle1)).*ones(1,length(taxis)),'w--');hold on;
% plot(taxis,(0-min([VPangle1_upper,VPangle1_lower])-mean(VPangle1)).*ones(1,length(taxis)),'w--');hold on;

% plot([50,50]./60,[-pi,pi],'r--');hold on;
% plot([250,250]./60,[-pi,pi],'r--');hold on;
% set(gca,'xtick',[50./60,(50+25*60)./60],'xticklabel',{'0','25'});hold on;
% xlim([0 (50+25*60)./60]);

% pbaspect([105*0.75 100 1]);
% xlim([25,55]); % 26-11:26+14

set(gca,'ytick',[-pi,0,pi],'yticklabel',{'-\pi','0','\pi'},'fontsize',14);
title(['\Deltakymo ' filename],'fontsize',14);
xlabel('Time (min)','fontsize',16);
ylabel('Polar angle \phi','fontsize',16);


% ylim([((0-min([VPangle1_upper,VPangle1_lower])-mean(VPangle1))+(0+max([VPangle1_upper,VPangle1_lower])-mean(VPangle1)))/2-pi/2,...
%     ((0-min([VPangle1_upper,VPangle1_lower])-mean(VPangle1))+(0+max([VPangle1_upper,VPangle1_lower])-mean(VPangle1)))/2+pi/2]);

% set(gca,'ytick',((0-min([VPangle1_upper,VPangle1_lower])-mean(VPangle1))+(0+max([VPangle1_upper,VPangle1_lower])-mean(VPangle1)))/2-pi/2+[-pi/2,0,pi/2]);

colormap(colorcet('L01'));
% daspect([7,1,1]);

% caxis([-3 9]);

% plot([taxis(t_on) taxis(t_on)],[-pi pi],'color','w','linewidth',1,'linestyle','--');hold on;
% plot([taxis(t_off) taxis(t_off)],[-pi pi],'color','w','linewidth',1,'linestyle','--');hold on;

plot([taxis(325) taxis(325)],[-pi pi],'color','w','linewidth',1,'linestyle','--');hold on;
plot([taxis(1000) taxis(1000)],[-pi pi],'color','w','linewidth',1,'linestyle','--');hold on;

rhos_interp(:,t_on);

col_cm = mean(bws{t_on}(:,2));row_cm = mean(bws{t_on}(:,2));

%%
isinROI1 = inpolygon(row_cm+rhos_interp(:,t_on).*sin(phis_interpcor_minuspitopi'+VPangle*pi/180)./umperpix,...
    col_cm+rhos_interp(:,t_on).*cos(phis_interpcor_minuspitopi'+VPangle*pi/180)./umperpix,...
    [ROI1(2,1),ROI1(2,2),ROI1(2,2),ROI1(2,1)],[ROI1(1,1),ROI1(1,1),ROI1(1,2),ROI1(1,2)]);

isinROI2 = inpolygon(row_cm+rhos_interp(:,t_on).*sin(phis_interpcor_minuspitopi'+VPangle*pi/180)./umperpix,...
    col_cm+rhos_interp(:,t_on).*cos(phis_interpcor_minuspitopi'+VPangle*pi/180)./umperpix,...
    [ROI2(2,1),ROI2(2,2),ROI2(2,2),ROI2(2,1)],[ROI2(1,1),ROI2(1,1),ROI2(1,2),ROI2(1,2)]);

%%
phis_ROI1 = find(isinROI1);phis_ROI2 = find(isinROI2);

plot([taxis(t_on) taxis(t_off)],[phis_interpcor_minuspitopi(phis_ROI1(1)) phis_interpcor_minuspitopi(phis_ROI1(1))],...
    'color','c','linewidth',1.5,'linestyle','-');hold on;
plot([taxis(t_on) taxis(t_off)],[phis_interpcor_minuspitopi(phis_ROI1(end)) phis_interpcor_minuspitopi(phis_ROI1(end))],...
    'color','c','linewidth',1.5,'linestyle','-');hold on;
plot([taxis(t_on) taxis(t_on)],[phis_interpcor_minuspitopi(phis_ROI1(1)) phis_interpcor_minuspitopi(phis_ROI1(end))],...
    'color','c','linewidth',1.5,'linestyle','-');hold on;
plot([taxis(t_off) taxis(t_off)],[phis_interpcor_minuspitopi(phis_ROI1(1)) phis_interpcor_minuspitopi(phis_ROI1(end))],...
    'color','c','linewidth',1.5,'linestyle','-');hold on;

plot([taxis(t_on) taxis(t_off)],[phis_interpcor_minuspitopi(phis_ROI2(1)) phis_interpcor_minuspitopi(phis_ROI2(1))],...
    'color','c','linewidth',1.5,'linestyle','-');hold on;
plot([taxis(t_on) taxis(t_off)],[phis_interpcor_minuspitopi(phis_ROI2(end)) phis_interpcor_minuspitopi(phis_ROI2(end))],...
    'color','c','linewidth',1.5,'linestyle','-');hold on;
plot([taxis(t_on) taxis(t_on)],[phis_interpcor_minuspitopi(phis_ROI2(1)) phis_interpcor_minuspitopi(phis_ROI2(end))],...
    'color','c','linewidth',1.5,'linestyle','-');hold on;
plot([taxis(t_off) taxis(t_off)],[phis_interpcor_minuspitopi(phis_ROI2(1)) phis_interpcor_minuspitopi(phis_ROI2(end))],...
    'color','c','linewidth',1.5,'linestyle','-');hold on;

pbaspect([1,1,1]);

caxis([-50 120]);

fig_current = gcf; fig_current.Renderer = 'painters';
print(fig_current,['kymo_' filename '_ROI1_2_cetL01'],'-dpdf');


%%

% %% VP1, curv local 
% 
% VPangle = mean(VPangle1*180/pi);% in degree
% shiftVP = findequal(phis_interpcor,VPangle*pi/180,1);
% 
% kymoss_interp = [curvs_interp;curvs_interp;curvs_interp];
% kymos_interp_shiftVP = kymoss_interp(N_phiinterp+shiftVP-N_phiinterp/2:N_phiinterp+shiftVP+N_phiinterp/2,:);
% 
% kymos_interp_shiftVP0 = mean(kymos_interp_shiftVP(:,1:50),2)*ones(1,t_total);
% 
% figure;
% 
% imagesc(taxis,phis_interpcor,kymos_interp_shiftVP-kymos_interp_shiftVP0);colorbar;hold on;
% plot(taxis,(0+max([VPangle1_upper,VPangle1_lower])-mean(VPangle1)).*ones(1,length(taxis)),'w--');hold on;
% plot(taxis,(0-min([VPangle1_upper,VPangle1_lower])-mean(VPangle1)).*ones(1,length(taxis)),'w--');hold on;
% 
% plot([50,50]./60,[-pi,pi],'r--');hold on;
% plot([250,250]./60,[-pi,pi],'r--');hold on;
% set(gca,'xtick',[50./60,(50+25*60)./60],'xticklabel',{'0','25'});hold on;
% xlim([0 (50+25*60)./60]);
% 
% % set(gca,'ytick',[-pi,0,pi],'yticklabel',{'-\pi','0','\pi'},'fontsize',14);
% title(['\Deltacurv ' filename],'fontsize',14);
% xlabel('Time (min)','fontsize',16);
% ylabel('Polar angle \phi','fontsize',16);
% 
% 
% ylim([((0-min([VPangle1_upper,VPangle1_lower])-mean(VPangle1))+(0+max([VPangle1_upper,VPangle1_lower])-mean(VPangle1)))/2-pi/2,...
%     ((0-min([VPangle1_upper,VPangle1_lower])-mean(VPangle1))+(0+max([VPangle1_upper,VPangle1_lower])-mean(VPangle1)))/2+pi/2]);
% 
% set(gca,'ytick',((0-min([VPangle1_upper,VPangle1_lower])-mean(VPangle1))+(0+max([VPangle1_upper,VPangle1_lower])-mean(VPangle1)))/2-pi/2+[-pi/2,0,pi/2]);
% 
% colormap(colorcet('D09'));
% daspect([7,1,1]);
% 
% caxis([-0.001 0.001]);
% 
% % fig_current = gcf; fig_current.Renderer = 'painters';
% % print(fig_current,['curv_' filename 'ROI1_cetD09'],'-dpdf');
% 
% 
% %% VP2, kymo local
% 
% VPangle = mean(VPangle2*180/pi);% in degree
% shiftVP = findequal(phis_interpcor,VPangle*pi/180,1);
% 
% kymoss_interp = [kymos_interp;kymos_interp;kymos_interp];
% kymos_interp_shiftVP = kymoss_interp(N_phiinterp+shiftVP-N_phiinterp/2:N_phiinterp+shiftVP+N_phiinterp/2,:);
% 
% % figure;
% % imagesc(taxis,phis_interpcor,kymos_interp_shiftVP);colorbar;pbaspect([2,1,1]);hold on;
% % 
% % plot(taxis,(0+max([VPangle1_upper,VPangle1_lower])-mean(VPangle1)).*ones(1,length(taxis)),'w--');hold on;
% % plot(taxis,(0-min([VPangle1_upper,VPangle1_lower])-mean(VPangle1)).*ones(1,length(taxis)),'w--');hold on;
% % 
% % set(gca,'ytick',[-pi,0,pi],'yticklabel',{'-\pi','0','\pi'},'fontsize',14);
% % title(['kymo ' filename],'fontsize',14);
% % xlabel('Time (min)','fontsize',16); 
% % ylabel('Polar angle \phi','fontsize',16);
% 
% 
% 
% kymos_interp_shiftVP0 = mean(kymos_interp_shiftVP(:,1151:1200),2)*ones(1,t_total);
% 
% % kymos_interp_shiftVP0 = kymos_interp_shiftVP; %[-(kymos_interp_shiftVP(:,74)-kymos_interp_shiftVP(:,73))*ones(1,73) zeros(size(kymos_interp_shiftVP,1),size(kymos_interp_shiftVP,2)-73)];
% % for i = 1:1:size(kymos_interp_shiftVP0,2)
% %     kymos_interp_shiftVP(:,i) = kymos_interp_shiftVP0(:,i)./mean(kymos_interp_shiftVP0(:,i));
% % end
% 
% 
% figure;
% 
% imagesc(taxis,phis_interpcor,kymos_interp_shiftVP-kymos_interp_shiftVP0);colorbar;hold on; % this is kymo
% plot(taxis,(0+max([VPangle2_upper,VPangle2_lower])-mean(VPangle2)).*ones(1,length(taxis)),'w--');hold on;
% plot(taxis,(0+min([VPangle2_upper,VPangle2_lower])-mean(VPangle2)).*ones(1,length(taxis)),'w--');hold on;
% 
% plot([1200,1200]./60,[-pi,pi],'r--');hold on;
% plot([1400,1400]./60,[-pi,pi],'r--');hold on;
% set(gca,'xtick',[1200./60,(1200+25*60)./60],'xticklabel',{'0','25'});hold on;
% xlim([1150/60 (1200+25*60)./60]);
% 
% % set(gca,'ytick',[-pi,0,pi],'yticklabel',{'-\pi','0','\pi'},'fontsize',14);
% title(['\Deltakymo ' filename],'fontsize',14);
% xlabel('Time (min)','fontsize',16);
% ylabel('Polar angle \phi','fontsize',16);
% 
% 
% ylim([((0+min([VPangle2_upper,VPangle2_lower])-mean(VPangle2))+(0+max([VPangle2_upper,VPangle2_lower])-mean(VPangle2)))/2-pi/2,...
%     ((0+min([VPangle2_upper,VPangle2_lower])-mean(VPangle2))+(0+max([VPangle2_upper,VPangle2_lower])-mean(VPangle2)))/2+pi/2]);
% 
% set(gca,'ytick',((0+min([VPangle2_upper,VPangle2_lower])-mean(VPangle2))+(0+max([VPangle2_upper,VPangle2_lower])-mean(VPangle2)))/2-pi/2+[-pi/2,0,pi/2]);
% 
% colormap(colorcet('L01'));
% daspect([7,1,1]);
% 
% caxis([-0.8 1.2]);
% 
% % fig_current = gcf; fig_current.Renderer = 'painters';
% % print(fig_current,['kymo_' filename 'ROI2_cetL01'],'-dpdf');
% 
% %% VP2, curv local
% 
% VPangle = mean(VPangle2*180/pi);% in degree
% shiftVP = findequal(phis_interpcor,VPangle*pi/180,1);
% 
% kymoss_interp = [curvs_interp;curvs_interp;curvs_interp];
% kymos_interp_shiftVP = kymoss_interp(N_phiinterp+shiftVP-N_phiinterp/2:N_phiinterp+shiftVP+N_phiinterp/2,:);
% 
% kymos_interp_shiftVP0 = mean(kymos_interp_shiftVP(:,1151:1200),2)*ones(1,t_total);
% 
% figure;
% 
% imagesc(taxis,phis_interpcor,kymos_interp_shiftVP-kymos_interp_shiftVP0);colorbar;hold on; % this is kymo
% plot(taxis,(0+max([VPangle2_upper,VPangle2_lower])-mean(VPangle2)).*ones(1,length(taxis)),'w--');hold on;
% plot(taxis,(0+min([VPangle2_upper,VPangle2_lower])-mean(VPangle2)).*ones(1,length(taxis)),'w--');hold on;
% 
% plot([1200,1200]./60,[-pi,pi],'r--');hold on;
% plot([1400,1400]./60,[-pi,pi],'r--');hold on;
% set(gca,'xtick',[1200./60,(1200+25*60)./60],'xticklabel',{'0','25'});hold on;
% xlim([1150/60 (1200+25*60)./60]);
% 
% % set(gca,'ytick',[-pi,0,pi],'yticklabel',{'-\pi','0','\pi'},'fontsize',14);
% title(['\Deltacurv ' filename],'fontsize',14);
% xlabel('Time (min)','fontsize',16);
% ylabel('Polar angle \phi','fontsize',16);
% 
% 
% ylim([((0+min([VPangle2_upper,VPangle2_lower])-mean(VPangle2))+(0+max([VPangle2_upper,VPangle2_lower])-mean(VPangle2)))/2-pi/2,...
%     ((0+min([VPangle2_upper,VPangle2_lower])-mean(VPangle2))+(0+max([VPangle2_upper,VPangle2_lower])-mean(VPangle2)))/2+pi/2]);
% 
% set(gca,'ytick',((0+min([VPangle2_upper,VPangle2_lower])-mean(VPangle2))+(0+max([VPangle2_upper,VPangle2_lower])-mean(VPangle2)))/2-pi/2+[-pi/2,0,pi/2]);
% 
% colormap(colorcet('D09'));
% daspect([7,1,1]);
% 
% caxis([-0.001 0.001]);
% 
% % fig_current = gcf; fig_current.Renderer = 'painters';
% % print(fig_current,['curv_' filename 'ROI2_cetD09'],'-dpdf');





%% functions

function num = findequal(array,value,k)
nums = find(array-value==abs(array-value));
if k<=length(nums)
    num = nums(1:k);
else
    disp('demanded output exceeds located equal-value numbers in array!\n');
    num = nums;
end
end

function [xCoM,yCoM] = getxyCoM(bws_x,bws_y)
%% Calculate center of mass from boundary pixel coordinate list x and y
if length(bws_x)~=length(bws_y)
    disp('x and y bounday not same length!');
    xCoM = [];yCoM = [];
else
    t_total = length(bws_x);
    xCoM = zeros(1,t_total);
    yCoM = zeros(1,t_total);
    
    for t = 1:1:t_total
        weightArea = 0;
        xSum = 0;ySum = 0;
        x_N = [bws_x{t};bws_x{t}(1)];y_N = [bws_y{t};bws_y{t}(1)];
        for n = 1:1:length(bws_x{t})
            weightArea = weightArea + (x_N(n)*y_N(n+1)-x_N(n+1)*y_N(n));
            xSum = xSum + (x_N(n)+x_N(n+1))*(x_N(n)*y_N(n+1)-x_N(n+1)*y_N(n));
            ySum = ySum + (y_N(n)+y_N(n+1))*(x_N(n)*y_N(n+1)-x_N(n+1)*y_N(n));
        end
        weightArea = weightArea/2;
        xCoM(t) = xSum/(6*weightArea);
        yCoM(t) = ySum/(6*weightArea);
        
    end
end

end

function [rhos_itp,phis_itp,I_shift] = getpolarCoM(t_range,cols_itp,rows_itp,xCoM_itp,yCoM_itp,drawon) % convert to polar coordinate & in microns
% add output the shifted term

global Lmax umperpix

rhos_itp = cell(1,length(t_range));
phis_itp = cell(1,length(t_range));
I_shift = cell(1,length(t_range));

if drawon
    figure;
end

for i_range = 1:1:length(t_range)
    
    t = t_range(i_range); % frame
    
    disp(t);
    
    col_itp = cols_itp{t};row_itp = rows_itp{t};
    
    rho_t = sqrt((col_itp-xCoM_itp(t)).^2+(row_itp-yCoM_itp(t)).^2).*umperpix; % in microns
    phi_t = atan2(row_itp-yCoM_itp(t),col_itp-xCoM_itp(t));
    
    % figure;
    % plot([col_itp;col_itp(1)],[row_itp;row_itp(1)],'r-','linewidth',1.5);hold on;daspect([1,1,1]);
    % set(gca,'ydir','reverse');
    % scatter(xCoM_itp(t),yCoM_itp(t),30,'r','filled');hold on;
    % for i = 1:10:length(col_itp)
    %     plot([xCoM_itp(t),col_itp(i)],[yCoM_itp(t),row_itp(i)],'k--');hold on;
    %     text(col_itp(i),row_itp(i),[num2str(phi_t(i)*180/pi,'%.1f'),char(176)],'fontsize',14);hold on;
    % end
    
    [phi_t_sorted,I_t] = sort(phi_t,'ascend');
    rho_t_sorted = rho_t(I_t);
    
    rhos_itp{t} = rho_t_sorted;
    phis_itp{t} = phi_t_sorted;
    I_shift{t} = I_t;
    
    if drawon
        
        % figure;
        % plot(phi_t_sorted);hold on;
        
        plot((col_itp-xCoM_itp(t)).*umperpix,(row_itp-yCoM_itp(t)).*umperpix,'linewidth',1.5);hold on;
        plot(rho_t_sorted.*cos(phi_t_sorted),rho_t_sorted.*sin(phi_t_sorted),'r:','linewidth',2);hold on;
        daspect([1,1,1]);
        axis([-Lmax Lmax -Lmax Lmax].*umperpix);
        
        %%
        
        set(gca,'ydir','reverse');
        
        pause(0.02);
        
        if i_range<length(t_range)
            clf;
        end
    end
    
end

end
