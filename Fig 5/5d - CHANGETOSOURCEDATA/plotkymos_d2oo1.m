%% check bws from raw dic video

global spf umperpix Lmax size_col size_row

filename = 'd2oo1';

umperpix = 0.691;

graylim = [39 133];imArrlim = [1 109];

size_col = 512; size_row = 512;

Lmax = 150/umperpix;

t_totallist = {50,672};
    
spflist = {5,11.38};

suffixlist = {'_t1','_t2'};

dictifpath = {'/Volumes/One Touch/2022 Photoactivation - Analysis/Analysis_FastPhotoactivatableWaves_OptoEct2_20220411/d2oo1/d2oo1_lightoff_tseries1_spf5s.czi - DIC.tif',...
    '/Volumes/One Touch/2022 Photoactivation - Analysis/Analysis_FastPhotoactivatableWaves_OptoEct2_20220411/d2oo1/d2oo1_lighton_tseries2_spf12s_power15_ite3_globaltrack2on.czi - DIC.tif'};
geftifpath = {'/Volumes/One Touch/2022 Photoactivation - Analysis/Analysis_FastPhotoactivatableWaves_OptoEct2_20220411/d2oo1/d2oo1_lightoff_tseries1_spf5s.czi - C=0.tif',...
    '/Volumes/One Touch/2022 Photoactivation - Analysis/Analysis_FastPhotoactivatableWaves_OptoEct2_20220411/d2oo1/d2oo1_lighton_tseries2_spf12s_power15_ite3_globaltrack2on.czi - C=0.tif'};

%%

w = 8;sFac2 = 0.99; % contract extracted bws again to match optimal intensity kymograph extraction

kymos_t1_2 = zeros(3500,50+672);

load(['kymos_updated/kymos_' filename '_t1.mat'],'kymos');
kymos_t1_2(1:size(kymos,1),1:50) = kymos;

load(['kymos_updated/kymos_' filename '_t2.mat'],'kymos');
kymos_t1_2(1:size(kymos,1),51:672+50) = kymos;
    

kymos_t1_2 = kymos_t1_2(:,1:(672+50));

figure;
imagesc(kymos_t1_2);colorbar;caxis(imArrlim);hold on;pbaspect([2,1,1]);

bws_x_t1_2 = cell(1,672+50);bws_y_t1_2 = cell(1,672+50);
load(['bws_updated/bws_' filename suffixlist{1} '.mat'],'bws_x','bws_y');
for t = 1:1:t_totallist{1}
    bws_x_t1_2{t} = bws_x{t};
    bws_y_t1_2{t} = bws_y{t};
end
load(['bws_updated/bws_' filename suffixlist{2} '.mat'],'bws_x','bws_y');
for t = (1+t_totallist{1}):1:(672+50)
    bws_x_t1_2{t} = bws_x{t-t_totallist{1}};
    bws_y_t1_2{t} = bws_y{t-t_totallist{1}};
end

%%

t_total = 672+50;

taxis = (0:1:(t_totallist{1}-1)).*spflist{1}/60;
for k_file = 2:1:length(suffixlist)
    
    taxis = [taxis taxis(end)+(1:1:(t_totallist{k_file})).*spflist{k_file}/60]; % in min
    
end

taxis = taxis(1:t_total);

%% change to polar coordinate and visualize

[xCoM_raw,yCoM_raw] = getxyCoM(bws_x_t1_2,bws_y_t1_2); % in pixels

t = 1:1:t_total;

drawon = 0;
[rhos_raw,phis_raw,I_shift] = getpolarCoM(t,bws_x_t1_2,bws_y_t1_2,xCoM_raw,yCoM_raw,drawon); % convert to polar coordinate & in microns

kymos_t1_2_shifted = nan(size(kymos));
for t = 1:1:t_total
    N_bd = length(bws_x_t1_2{t});
    kymo = kymos_t1_2(1:N_bd,t);
    
    kymos_t1_2_shifted(1:N_bd,t) = kymo(I_shift{t});
end

%% compare interpolate using interp1 (interp twice on 180 degree apart breakings on the circle and combine)

N_phiinterp = 36000; % 0.1 degree
phis_interpcor_minuspitopi = linspace(-pi,pi,N_phiinterp+1);
phis_interpcor_0to2pi = linspace(0,2*pi,N_phiinterp+1);

kymos_interp_minuspitopi = zeros(N_phiinterp+1,t_total); 
kymos_interp_0to2pi = zeros(N_phiinterp+1,t_total);

for t = 1:1:t_total
    
    N_bd = length(bws_x_t1_2{t});
    
    phit = phis_raw{t};
    [uniq_phit,I_phi] = unique(phit);
    kymo = kymos_t1_2_shifted(1:N_bd,t);
    kymo_interp = interp1(uniq_phit,kymo(I_phi),phis_interpcor_minuspitopi);
    kymos_interp_minuspitopi(:,t) = kymo_interp;
    
    phit2 = phis_raw{t};
    phit2(phit2<0) = phit2(phit2<0)+2*pi;
    [uniq_phit,I_phi] = unique(phit2);
    kymo = kymos_t1_2_shifted(1:N_bd,t);
    kymo_interp = interp1(uniq_phit,kymo(I_phi),phis_interpcor_0to2pi);
    kymos_interp_0to2pi(:,t) = kymo_interp;
    
end

kymos_interp = nan(N_phiinterp+1,t_total);% on -\pi to \pi
phis_interpcor= linspace(-pi,pi,N_phiinterp+1);
kymos_interp1 = kymos_interp_minuspitopi;
kymos_interp2 = kymos_interp_0to2pi;kymos_interp2 = [kymos_interp2(N_phiinterp/2+1:N_phiinterp+1,:);kymos_interp2(1:N_phiinterp/2,:)];

mask_kymo1nan = (isnan(kymos_interp1))&(~isnan(kymos_interp2));
mask_kymo2nan = (~isnan(kymos_interp1))&(isnan(kymos_interp2));
mask_nonan = (~isnan(kymos_interp1))&(~isnan(kymos_interp2));
kymos_interp(mask_kymo1nan) = kymos_interp2(mask_kymo1nan);
kymos_interp(mask_kymo2nan) = kymos_interp1(mask_kymo2nan);
kymos_interp(mask_nonan) = (kymos_interp1(mask_nonan)+kymos_interp1(mask_nonan)).*0.5;


figure;
imagesc(kymos_interp_minuspitopi);colorbar;pbaspect([2,1,1]);hold on;
title('kymo interp -\pi to \pi','fontsize',16);
figure;
imagesc(kymos_interp_0to2pi);colorbar;pbaspect([2,1,1]);hold on;
title('kymo interp 0 to 2\pi','fontsize',16);
figure;
imagesc(kymos_interp);colorbar;pbaspect([2,1,1]);hold on;
title('kymo interp combined (-\pi to \pi)','fontsize',16);

%% periodic and shift center

VPangle = -48;% in degree, -42 to start with
shiftVP = findequal(phis_interpcor,VPangle*pi/180,1);

kymoss_interp = [kymos_interp;kymos_interp;kymos_interp];
kymos_interp_shiftVP = kymoss_interp(N_phiinterp+shiftVP-N_phiinterp/2:N_phiinterp+shiftVP+N_phiinterp/2,:);

figure;
imagesc(taxis,phis_interpcor,kymos_interp_shiftVP);colorbar;pbaspect([2,1,1]);hold on;
caxis([40 110]);
set(gca,'ytick',[-pi,0,pi],'yticklabel',{'-\pi','0','\pi'},'fontsize',14);
title(['kymo ' filename],'fontsize',14);
xlabel('Time (min)','fontsize',16);
ylabel('Polar angle \phi','fontsize',16);
daspect([7,1,1]);

%%

% kymos_interp_shiftVP0 =
% mean(kymos_interp_shiftVP(:,1:50),2)*ones(1,t_total); % to start with

kymos_interp_shiftVP0 = kymos_interp_shiftVP; %[-(kymos_interp_shiftVP(:,74)-kymos_interp_shiftVP(:,73))*ones(1,73) zeros(size(kymos_interp_shiftVP,1),size(kymos_interp_shiftVP,2)-73)];

for i = 1:1:size(kymos_interp_shiftVP0,2)
    kymos_interp_shiftVP(:,i) = kymos_interp_shiftVP0(:,i)./mean(kymos_interp_shiftVP0(:,i));
end


figure;
imagesc(taxis(51:end),phis_interpcor,kymos_interp_shiftVP(:,51:end));colorbar;hold on;

plot([taxis(242+50) taxis(242+50)],[-pi pi],'color','w','linewidth',1,'linestyle','--');hold on;
plot([taxis(253+50) taxis(253+50)],[-pi pi],'color','w','linewidth',1,'linestyle','--');hold on;
plot([taxis(264+50) taxis(264+50)],[-pi pi],'color','w','linewidth',1,'linestyle','--');hold on;

% caxis([-30 50]);
caxis([0.5 1.5]);
set(gca,'ytick',[-pi,0,pi],'yticklabel',{'-\pi','0','\pi'},'fontsize',14);
title(['\Deltakymo ' filename],'fontsize',14);
xlabel('Time (min)','fontsize',16);
ylabel('Polar angle \phi','fontsize',16);

colormap(colorcet('L01'));
% daspect([7,1,1]);

pbaspect([26.25*(taxis(end)-taxis(51))/10 100 1]);

% xlim([50-12,50+13]); % 50-12:50+13

% set(gca,'xtick',[38,43,48,53,58,63,68]);

% xlim([27.5,27.5]+[0,30]);

fig_current = gcf; fig_current.Renderer = 'painters';
print(fig_current,['normkymo_' filename '_cetL01_20230218'],'-dpdf');

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