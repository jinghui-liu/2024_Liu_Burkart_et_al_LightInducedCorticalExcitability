%% check bws from raw dic video

clear all

global spf umperpix Lmax size_col size_row

filename = 'optoect220220307';

umperpix = 266.74/512;

rholim = [0 200];

size_col = 512; size_row = 512;

Lmax = 150/umperpix;

t_total = 1131;

spf = 8;

taxis = (0:1:(t_total-1)).*spf/60;

% dictifpath = '/Volumes/One Touch/2022 Photoactivation - Analysis/Analysis_WildtypeSCW/d7oo1-14wt/d7oo11wt_dic8bit.tif';
geftifpath = '/Volumes/One Touch/2022 Photoactivation - Analysis/2 Analysis_ReversibleContractions/optoect220220307/d2oo2_optoect220220307.tif';

% load(['bws_updated/bws_' filename '.mat'],'bws_x','bws_y','rad','thresh','smfac','sFac');

load('save_kymocurvextract_optoect220220307_20220911 17_11/boundary.mat','bws');

bws_x = cell(1,t_total);bws_y = cell(1,t_total);
for t = 1:1:t_total
    bws_x{t} = bws{t}(:,2);
    bws_y{t} = bws{t}(:,1);
end

%% change to polar coordinate and visualize

% load(['curvs_updated/curvs_' filename '.mat'],'curvs','arcpoints','Ninterpcurv');

load('save_kymocurvextract_optoect220220307_20220911 17_11/curvs.mat','curvs');

% figure;imagesc(curvs);colorbar;

[xCoM_raw,yCoM_raw] = getxyCoM(bws_x,bws_y); % in pixels

t = 1:1:t_total;

drawon = 0;
[rhos_raw,phis_raw,I_shift] = getpolarCoM(t,bws_x,bws_y,xCoM_raw,yCoM_raw,drawon); % convert to polar coordinate & in microns

curvs_shifted = nan(size(curvs));
for t = 1:1:t_total
    N_bd = length(bws_x{t});
    curv = curvs(1:N_bd,t);
    
    curvs_shifted(1:N_bd,t) = curv(I_shift{t});
end



%% compare interpolate using interp1 (interp twice on 180 degree apart breakings on the circle and combine)

N_phiinterp = 3600; % 0.01 degree
phis_interpcor_minuspitopi = linspace(-pi,pi,N_phiinterp+1);
phis_interpcor_0to2pi = linspace(0,2*pi,N_phiinterp+1);

curvs_interp_minuspitopi = zeros(N_phiinterp+1,t_total); 
curvs_interp_0to2pi = zeros(N_phiinterp+1,t_total);

rhos_interp_minuspitopi = zeros(N_phiinterp+1,t_total); 
rhos_interp_0to2pi = zeros(N_phiinterp+1,t_total);

for t = 1:1:t_total
    
    N_bd = length(bws_x{t});
    
    phit = phis_raw{t};
    [uniq_phit,I_phi] = unique(phit);
    curv = curvs_shifted(1:N_bd,t);
    curv_interp = interp1(uniq_phit,curv(I_phi),phis_interpcor_minuspitopi);
    curvs_interp_minuspitopi(:,t) = curv_interp;
    
    rhot = rhos_raw{t}; % in microns
    rhot_interp = interp1(uniq_phit,rhot(I_phi),phis_interpcor_minuspitopi);
    rhos_interp_minuspitopi(:,t) = rhot_interp;
    
    phit2 = phis_raw{t};
    phit2(phit2<0) = phit2(phit2<0)+2*pi;
    [uniq_phit,I_phi] = unique(phit2);
    curv = curvs_shifted(1:N_bd,t);
    curv_interp = interp1(uniq_phit,curv(I_phi),phis_interpcor_0to2pi);
    curvs_interp_0to2pi(:,t) = curv_interp;
    
    rhot2 = rhos_raw{t}; % in microns
    rhot_interp = interp1(uniq_phit,rhot(I_phi),phis_interpcor_0to2pi);
    rhos_interp_0to2pi(:,t) = rhot_interp;
    
end

phis_interpcor= linspace(-pi,pi,N_phiinterp+1);

curvs_interp = nan(N_phiinterp+1,t_total);% on -\pi to \pi
curvs_interp1 = curvs_interp_minuspitopi;
curvs_interp2 = curvs_interp_0to2pi;curvs_interp2 = [curvs_interp2(N_phiinterp/2+1:N_phiinterp+1,:);curvs_interp2(1:N_phiinterp/2,:)];

mask_curv1nan = (isnan(curvs_interp1))&(~isnan(curvs_interp2));
mask_curv2nan = (~isnan(curvs_interp1))&(isnan(curvs_interp2));
mask_nonan = (~isnan(curvs_interp1))&(~isnan(curvs_interp2));
curvs_interp(mask_curv1nan) = curvs_interp2(mask_curv1nan);
curvs_interp(mask_curv2nan) = curvs_interp1(mask_curv2nan);
curvs_interp(mask_nonan) = (curvs_interp1(mask_nonan)+curvs_interp1(mask_nonan)).*0.5;

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
% imagesc(curvs_interp_minuspitopi);colorbar;pbaspect([2,1,1]);hold on;
% title('curv interp -\pi to \pi','fontsize',16);
% figure;
% imagesc(curvs_interp_0to2pi);colorbar;pbaspect([2,1,1]);hold on;
% title('curv interp 0 to 2\pi','fontsize',16);
% figure;
% imagesc(curvs_interp);colorbar;pbaspect([2,1,1]);hold on;
% title('curv interp combined (-\pi to \pi)','fontsize',16);

%% periodic and shift center

VPangle = 108-180;% in degree
shiftVP = findequal(phis_interpcor,VPangle*pi/180,1);

curvss_interp = [curvs_interp;curvs_interp;curvs_interp];
curvs_interp_shiftVP = curvss_interp(N_phiinterp+shiftVP-N_phiinterp/2:N_phiinterp+shiftVP+N_phiinterp/2,:);

figure;
imagesc(taxis,phis_interpcor,curvs_interp_shiftVP);colorbar;pbaspect([2,1,1]);hold on; % in inverse pixel

colormap(colorcet('D09'));

caxis([0.003 0.006]);

set(gca,'ytick',[-pi,0,pi],'yticklabel',{'-\pi','0','\pi'},'fontsize',14);
title(['curv ' filename],'fontsize',14);
xlabel('Time (min)','fontsize',16);
ylabel('Polar angle \phi','fontsize',16);

%%

curvs_interp_shiftVP0 = mean(curvs_interp_shiftVP(:,1:1),2)*ones(1,t_total); % 61:1:61

figure;
imagesc(taxis,phis_interpcor,(curvs_interp_shiftVP-curvs_interp_shiftVP0)./umperpix);cbh = colorbar;hold on; % in inverse micron

colormap(colorcet('D09'));
caxis([-0.002 0.002]);
ylabel(cbh,'Curvature (\mum^{-1})','fontsize',16);

set(gca,'ytick',[-pi,0,pi],'yticklabel',{'-\pi','0','\pi'},'fontsize',14);
title(['\Deltacurv ' filename],'fontsize',14);
xlabel('Time (min)','fontsize',16);
ylabel('Polar angle \phi','fontsize',16);

% daspect([7,1,1]);
% pbaspect([105*0.75 100 1]);
% xlim([25,55]); % 26-11:26+14

% set(gca,'xtick',20:10:60);

% fig_current = gcf; fig_current.Renderer = 'painters';
% print(fig_current,['deltacurv_' filename '_cetD09_highthres'],'-dpdf');

curv_plotted = (curvs_interp_shiftVP-curvs_interp_shiftVP0)./umperpix;

% plot([taxis(73) taxis(73)],[-pi pi],'color','k','linewidth',1.5,'linestyle','--');hold on;
% plot([taxis(89) taxis(89)],[-pi pi],'color','k','linewidth',1.5,'linestyle','--');hold on;
% plot([taxis(105) taxis(105)],[-pi pi],'color','k','linewidth',1.5,'linestyle','--');hold on;

% %% Gaussian peak fitting
% 
% phimins = cell(1,t_total);
% taxismins = cell(1,t_total);
% 
% % savedir = 'd7oo11wt_plotcurvs_fitpks';
% % mkdir(savedir);
% 
% % figure;
% 
% for t = 1:1:t_total
%     
%     curv = smooth(curv_plotted(:,t),50);phi = phis_interpcor;
%     
%     locs = find(curv<0);
%     
%     if isempty(locs)
%         
%         phimins{t} = [];
%         taxismins{t} = [];
%         
%     else
%         
%         loc_seg = cell(1,1);
%         i_seg = 1;
%         loc_seg{1}(1) = locs(1);
%         for i = 1:1:length(locs)-1
%             if locs(i+1)-locs(i)>1
%                 loc_seg{i_seg}(2) = locs(i);
%                 loc_seg{i_seg+1}(1) = locs(i+1);
%                 i_seg = i_seg+1;
%             end
%         end
%         loc_seg{i_seg}(2) =locs(end);
%         
%         loc_seg2 = cell(1,1);
%         i_seg2 = 1;
%         for i_seg = 1:1:length(loc_seg)
%             if loc_seg{i_seg}(2)-loc_seg{i_seg}(1)>100 % customized cutoff
%                 loc_seg2{i_seg2} = loc_seg{i_seg};
%                 i_seg2 = i_seg2+1;
%             end
%         end
%         
%         if isempty(loc_seg2{1})
%             
%             phimins{t} = [];
%             taxismins{t} = [];
%         else
%             
% %             plot(phi,curv,'k-','linewidth',0.5);hold on;
% %             
% %             plot([-pi,pi],[0,0],'k--','linewidth',1.5);hold on;
%             
%             for i_seg2 = 1:1:length(loc_seg2)
%                 
%                 xdata = reshape(phi(loc_seg2{i_seg2}(1):loc_seg2{i_seg2}(2)),1,[]);
% %                 
%                 ydata = reshape(curv(loc_seg2{i_seg2}(1):loc_seg2{i_seg2}(2)),1,[]);
% %                 
% %                 plot(xdata,ydata,'b--','linewidth',1);hold on;
% %                 
% %                 fun = @(x,xdata)x(1)*exp(-((xdata-x(2))./x(3)).^2);
% %                 x0 = [1,0,1];
% %                 
% %                 x = lsqcurvefit(fun,x0,xdata,ydata);
% %                 
% %                 plot(xdata,fun(x,xdata),'r-','linewidth',1.5);hold on;
%                 
% %                 phimins{t}(i_seg2) = x(2);
% 
% [ydata_sorted,I] = sort(ydata,'ascend');
% 
% %                 phimins{t}(i_seg2) = mean(xdata(I(1:10)));
% 
% pos = loc_seg2{i_seg2}(1):loc_seg2{i_seg2}(2);
% 
%                 phimins{t}(i_seg2) = round(mean(pos(I(1:5))));
%                 taxismins{t}(i_seg2) = t;
%                 
%             end
%             
% %             xlim([-pi pi]);
% %             pbaspect([1.25 1 1]);
% %             
% %             xlabel('Polar angle \theta','fontsize',16);
% %             ylabel('\Delta curvature (\mum^{-2})','fontsize',16);
% %             
% %             title(['frame ' num2str(t) ', ' num2str(round((t-1)*spf/60)) ' min'],'fontsize',14);
% %             
% %             %%
% %             
% %             saveas(gca,[savedir '/t' num2str(t) '.png']);
% %             clf;
%         end
%     end
%     
% end
% 
% 
% figure;
% imagesc(taxis,phis_interpcor,(curvs_interp_shiftVP-curvs_interp_shiftVP0)./umperpix);cbh = colorbar;hold on; % in inverse micron
% 
% colormap(colorcet('D09'));
% caxis([-0.002 0.002]);
% ylabel(cbh,'Curvature (\mum^{-1})','fontsize',16);
% 
% set(gca,'ytick',[-pi,0,pi],'yticklabel',{'-\pi','0','\pi'},'fontsize',14);
% title(['\Deltacurv ' filename],'fontsize',14);
% xlabel('Time (min)','fontsize',16);
% ylabel('Polar angle \phi','fontsize',16);
% 
% xlim([25,55]);
% pbaspect([105*0.75 100 1]);
% 
% plot([taxis(73) taxis(73)],[-pi pi],'color','k','linewidth',1,'linestyle','--');hold on;
% plot([taxis(89) taxis(89)],[-pi pi],'color','k','linewidth',1,'linestyle','--');hold on;
% plot([taxis(105) taxis(105)],[-pi pi],'color','k','linewidth',1,'linestyle','--');hold on;
% 
% 
% % for t = 75:1:107
% %     for i = 1:1:length(phimins{t})
% %         scatter(taxis(t),phi(phimins{t}(i)),100,'w');hold on;
% %     end
% % end
% 
% 
% kymo_plotted = kymos_interp_shiftVP-kymos_interp_shiftVP0;
% 
% 
% A = [phimins{:}];
% T = [taxismins{:}];
% 
% A2 = [];
% T2 = [];
% 
% i2 = 0;
% for i = 1:1:length(A)
%     if (T(i)==103)&&(abs(A(i)-1800)>200)
% %         scatter(taxis(T(i)),phi(A(i)),100,'c','filled');hold on;
%         
%         i2 = i2+1;
%         T2(i2) = T(i);
%         A2(i2) = A(i);
%     else
%         if ((T(i)>=75)&&(T(i)<=103))&&(kymo_plotted(A(i),T(i))>1)
% %             scatter(taxis(T(i)),phi(A(i)),100,'c','filled');hold on;
%             
%             i2 = i2+1;
%             T2(i2) = T(i);
%             A2(i2) = A(i);
%         end
%     end
% end
% 
% [sorted_A2,I] = sort(A2,'ascend');
% sorted_T2 = T2(I);
% 
% plot(smooth(taxis(sorted_T2),5),smooth(phi(sorted_A2),5),'color','c','linewidth',1.5);hold on;
% 
% % Aq = linspace(min(A2),max(A2),21);
% % [~,I2] = unique(sorted_A2);
% % Tq = interp1(sorted_A2(I2),sorted_T2(I2),Aq);
% % plot(smooth(taxis(round(Tq)),5),smooth(phi(round(Aq)),5),'color','c','linewidth',1.5);hold on;
% 
% fig_current = gcf; fig_current.Renderer = 'painters';
% print(fig_current,['deltacurv_mincurvannotated_' '_cetD09_highthres'],'-dpdf');

%% Plot curv rings overlayed to outter oocyte

% imArr = readfromtiff(filename,1131,512,512);
% save(['imArr_' filename,'.mat'],'imArr','-v7.3');

load(['imArr_' filename,'.mat'],'imArr');

%%

savedir = [filename '_plotcurvs'];
mkdir(savedir);

lag_frame = 1;

[mycurvmap, ~,~] = colorcet('D09'); % visualize curvature contour around boundary

deltacurvs_interp = (curvs_interp-mean(curvs_interp(:,1:1),2)*ones(1,t_total))./umperpix;

dilate_fac = 1.125;

mincurv = -0.002; maxcurv = 0.002;

figure;

for t = 661:1:661 % 1:lag_frame:size(imArr,3) % 100,325,1000 
    
    disp(t);
    
    I = imArr(:,:,t);
    
    col_bd = bws{t}(:,2);row_bd = bws{t}(:,1);
    
    col_cm = mean(col_bd);row_cm = mean(row_bd);
    
    dilated_col_bd = col_cm+(col_bd-col_cm).*dilate_fac; % approx
    dilated_row_bd = row_cm+(row_bd-row_cm).*dilate_fac; % approx
    
    rhot_interp = rhos_interp(:,t)./umperpix; % transform from microns to pixels
    
    dilated_col_bd_interp = col_cm+dilate_fac.*rhot_interp.*cos(phis_interpcor_minuspitopi'); %%
    dilated_col_bd_interp_sm = smooth(dilated_col_bd_interp,20);
    
    dilated_row_bd_interp = row_cm+dilate_fac.*rhot_interp.*sin(phis_interpcor_minuspitopi'); %%
    dilated_row_bd_interp_sm = smooth(dilated_row_bd_interp,20);
    
%     if t==1
%         [X,Y] = meshgrid(1:1:size_col,1:1:size_row);
%         mask = inpolygon(X,Y,dilated_col_bd_interp_sm,dilated_row_bd_interp_sm); % get rid of fluorescence signals outside the circle
%     end
    %%
    
    imagesc(I);hold on;colormap(colorcet('L01'));caxis(rholim);
    
    curv = deltacurvs_interp(:,t);
    
    for i_curv = 1:1:length(dilated_col_bd_interp_sm)
        
        if curv(i_curv)<=mincurv
            colork = 1;
        else
            colork = 1+min(size(mycurvmap,1)-1,round(size(mycurvmap,1)*((curv(i_curv)-mincurv)/(maxcurv-mincurv))));
        end
        
        if i_curv<length(dilated_col_bd_interp_sm)
            %             plot(col_bd(i_curv:(i_curv+1)),row_bd(i_curv:(i_curv+1)),'color',curv_cm(colork,:),'linewidth',2);hold on;
%             plot(col_cm+(col_bd(i_curv:(i_curv+1))-col_cm).*dilate_fac,row_cm+(row_bd(i_curv:(i_curv+1))-row_cm).*dilate_fac,'color',mycurvmap(colork,:),'linewidth',2);hold on;

            plot([dilated_col_bd_interp_sm(i_curv),dilated_col_bd_interp_sm(i_curv+1)],[dilated_row_bd_interp_sm(i_curv),dilated_row_bd_interp_sm(i_curv+1)],'color',mycurvmap(colork,:),'linewidth',2);hold on;

        else
            %             plot(col_bd(i_curv):col_bd(1),row_bd(i_curv):row_bd(1),'color',curv_cm(colork,:),'linewidth',2);hold on;
%             plot(col_cm+(col_bd(i_curv):col_bd(1)-col_cm).*dilate_fac,row_cm+(row_bd(i_curv):row_bd(1)-row_cm).*dilate_fac,'color',mycurvmap(colork,:),'linewidth',2);hold on;

            plot([dilated_col_bd_interp_sm(i_curv),dilated_col_bd_interp_sm(1)],[dilated_row_bd_interp_sm(i_curv),dilated_row_bd_interp_sm(1)],'color',mycurvmap(colork,:),'linewidth',2);hold on;

        end
    end
    
    scatter(col_cm,row_cm,100,'r','filled');hold on;
    daspect([1,1,1]);
    
    set(gca,'xminortick','off','yminortick','off');
    set(gca,'xtick',[],'ytick',[]);
    
    title(['frame ' num2str(t) ', ' num2str(round((t-1)*spf/60)) ' min'],'fontsize',14);
    
    ylim([row_cm-150/umperpix row_cm+150/umperpix]);xlim([col_cm-150/umperpix col_cm+150/umperpix]);
    
    %%
%     saveas(gca,[savedir '/t' num2str(t) '.png']);
%     clf;

% fig_current = gcf; fig_current.Renderer = 'painters';
% print(fig_current,[savedir '/t' num2str(t) '_deltacurv_' '_cetD09_highthres'],'-dpdf');
    
end

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

