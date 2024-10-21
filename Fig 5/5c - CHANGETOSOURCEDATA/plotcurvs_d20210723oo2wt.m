global spf umperpix Lmax size_col size_row

filename = 'd20210723oo2wt';

umperpix = 0.521;

spf = 19;

graylim = [56 255];

size_col = 512; size_row = 512;

Lmax = 150/umperpix;

t_total = 303;

%%

taxis = (0:1:(t_total-1)).*spf;

load(['bws/bws_' filename '.mat'],'bws_x','bws_y','rad','thresh','smfac','sFac');

%% change to polar coordinate and visualize

load(['curvs/curvs_' filename '.mat'],'curvs','arcpoints','Ninterpcurv');

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

VPangle = 68;% in degree

shiftVP = findequal(phis_interpcor,VPangle*pi/180,1);

curvss_interp = [curvs_interp;curvs_interp;curvs_interp];
curvs_interp_shiftVP = curvss_interp(N_phiinterp+shiftVP-N_phiinterp/2:N_phiinterp+shiftVP+N_phiinterp/2,:);

figure;
imagesc(taxis,phis_interpcor,curvs_interp_shiftVP);colorbar;pbaspect([2,1,1]);hold on; % in inverse pixel

colormap(colorcet('D09'));

% caxis([0.003 0.006]);

set(gca,'ytick',[-pi,0,pi],'yticklabel',{'-\pi','0','\pi'},'fontsize',14);
title(['curv ' filename],'fontsize',14);
xlabel('Time (min)','fontsize',16);
ylabel('Polar angle \phi','fontsize',16);

%%

curvs_interp_shiftVP0 = mean(curvs_interp_shiftVP(:,1:6),2)*ones(1,length(taxis)); % 61:1:61

figure;
imagesc(taxis,phis_interpcor,(curvs_interp_shiftVP-curvs_interp_shiftVP0)./umperpix);cbh = colorbar;hold on; % in inverse micron

pbaspect([26.25*(taxis(end)/60)/10,100,1]);

% plot(taxis_plot,(0+max([VPangle1_upper,VPangle1_lower])-mean(VPangle1)).*ones(1,length(taxis_plot)),'w--','linewidth',1.5);hold on;
% plot(taxis_plot,(0+min([VPangle1_upper,VPangle1_lower])-mean(VPangle1)).*ones(1,length(taxis_plot)),'w--','linewidth',1.5);hold on;

colormap(colorcet('D09'));

clim([-0.003 0.003]);

ylabel(cbh,'Curvature (\mum^{-1})','fontsize',16);

set(gca,'ytick',[-pi,0,pi],'yticklabel',{'-\pi','0','\pi'},'fontsize',14);
title(['\Deltacurv ' filename],'fontsize',14);
xlabel('Time (s)','fontsize',16);
ylabel('Polar angle \phi','fontsize',16);


% mkdir('plot curvs');
plot([taxis(73) taxis(73)],[-pi pi],'color','k','linewidth',2,'linestyle','--');hold on; % frame = 73, 85, 97
plot([taxis(85) taxis(85)],[-pi pi],'color','k','linewidth',2,'linestyle','--');hold on;
plot([taxis(97) taxis(97)],[-pi pi],'color','k','linewidth',2,'linestyle','--');hold on;

curv_plotted = (curvs_interp_shiftVP-curvs_interp_shiftVP0)./umperpix;
% %
% fig_current = gcf; fig_current.Renderer = 'painters';
% print(fig_current,['plot curvs/deltacurv_' filename '_cetD09_highthres'],'-dpdf');

%%

% curvs_final = (curvs_interp_shiftVP(:,frame_axis_plot)-curvs_interp_shiftVP0)./umperpix;
%
% save(['curvs/' 'curvs_final_' filename '.mat'],'curvs_final','taxis_plot','phis_interpcor');

%% Gaussian peak fitting

phimins = cell(1,t_total);
taxismins = cell(1,t_total);

% savedir = [filename '_plotcurvs_fitpks'];
% mkdir(savedir);
%
% figure;

for t = 68:100

%     curv = smooth(curv_plotted(:,t),10);

    curv = curv_plotted(:,t);

    phi = phis_interpcor;

    locs = find(curv<0);

    if isempty(locs)

        phimins{t} = [];
        taxismins{t} = [];

    else

        loc_seg = cell(1,1);
        i_seg = 1;
        loc_seg{1}(1) = locs(1);
        for i = 1:1:length(locs)-1
            if locs(i+1)-locs(i)>1
                loc_seg{i_seg}(2) = locs(i);
                loc_seg{i_seg+1}(1) = locs(i+1);
                i_seg = i_seg+1;
            end
        end
        loc_seg{i_seg}(2) =locs(end);

        loc_seg2 = cell(1,1);
        i_seg2 = 1;
        for i_seg = 1:1:length(loc_seg)
            if loc_seg{i_seg}(2)-loc_seg{i_seg}(1)>5 % customized cutoff
                loc_seg2{i_seg2} = loc_seg{i_seg};
                i_seg2 = i_seg2+1;
            end
        end

        if isempty(loc_seg2{1})

            phimins{t} = [];
            taxismins{t} = [];
        else

            %             plot(phi,curv,'k-','linewidth',0.5);hold on;
            %
            %             plot([-pi,pi],[0,0],'k--','linewidth',1.5);hold on;

            for i_seg2 = 1:1:length(loc_seg2)

                xdata = reshape(phi(loc_seg2{i_seg2}(1):loc_seg2{i_seg2}(2)),1,[]);
                % %
                ydata = reshape(curv(loc_seg2{i_seg2}(1):loc_seg2{i_seg2}(2)),1,[]);
                % %
                %                 plot(xdata,ydata,'b--','linewidth',1);hold on;
                % %
                %                 fun = @(x,xdata)x(1)*exp(-((xdata-x(2))./x(3)).^2);
                %                 x0 = [1,0,1];
                %
                %                 x = lsqcurvefit(fun,x0,xdata,ydata);
                %
                %                 plot(xdata,fun(x,xdata),'r-','linewidth',1.5);hold on;
                %
                %                 phimins{t}(i_seg2) = x(2);

                [ydata_sorted,I] = sort(ydata,'ascend');

                %                 phimins{t}(i_seg2) = mean(xdata(I(1:10)));

                pos = loc_seg2{i_seg2}(1):loc_seg2{i_seg2}(2);

                phimins{t}(i_seg2) = round(mean(pos(I(1:5))));


                taxismins{t}(i_seg2) = t;

            end

            %             xlim([-pi pi]);
            %             pbaspect([1.25 1 1]);
            %
            %             xlabel('Polar angle \theta','fontsize',16);
            %             ylabel('\Delta curvature (\mum^{-2})','fontsize',16);
            %
            %             title(['frame ' num2str(t) ', ' num2str(round(taxis(t))) ' min'],'fontsize',14);

            %%

            %             saveas(gca,[savedir '/t' num2str(t) '.png']);
            %             clf;
        end
    end

end

%%

figure;
imagesc(taxis,phis_interpcor,curv_plotted);cbh = colorbar;hold on; % in inverse micron

colormap(colorcet('D09'));
clim([-0.003 0.003]);
ylabel(cbh,'Curvature (\mum^{-1})','fontsize',16);

plot([taxis(73) taxis(73)],[-pi pi],'color','k','linewidth',2,'linestyle','--');hold on; % frame = 73, 85, 97
plot([taxis(85) taxis(85)],[-pi pi],'color','k','linewidth',2,'linestyle','--');hold on;
plot([taxis(97) taxis(97)],[-pi pi],'color','k','linewidth',2,'linestyle','--');hold on;

set(gca,'ytick',[-pi,0,pi],'yticklabel',{'-\pi','0','\pi'},'fontsize',14);
title(['\Deltacurv ' filename],'fontsize',14);
xlabel('Time (min)','fontsize',16);
ylabel('Polar angle \phi','fontsize',16);

pbaspect([26.25*(taxis(end)/60)/10,100,1]);

% for t = 68:1:100
%     for i = 1:1:length(phimins{t})
%         scatter(taxis(t),phi(phimins{t}(i)),100,'w');hold on;
%     end
% end
% 
% A = [phimins{:}];
% T = [taxismins{:}];
% 
% for i = 1:1:length(A)
%     scatter(taxis(T(i)),phi(A(i)),100,'c','filled');hold on;
% end
% 
% % xlim([45 55]);
% 
% % manually delete noise points
% logical_A2_T2 = zeros(1,length(A));
% 
% t_ginput = zeros(1,100);phi_ginput = zeros(1,100);
% 
% for i = 1:1:100
% 
%     [t_ginput(i),phi_ginput(i)] = ginput(1);
% 
%     dist = sqrt((phi(A)-phi_ginput(i)).^2+(taxis(T)-t_ginput(i)).^2);
% 
%     pks = find(dist==min(dist));
% 
%     logical_A2_T2(pks) = 1;
% 
%     scatter(taxis(T(pks)),phi(A(pks)),120,'r');hold on;
% 
% end


load(['logical_A2_T2_' filename '.mat'],'logical_A2_T2');

T2 = T(logical(logical_A2_T2));
A2 = A(logical(logical_A2_T2));

[sorted_A2,I] = sort(A2,'ascend');
sorted_T2 = T2(I);
%
plot(smooth(taxis(sorted_T2),5),smooth(phi(sorted_A2),5),'color','c','linewidth',1.5);hold on;
%
% % Aq = linspace(min(A2),max(A2),21);
% % [~,I2] = unique(sorted_A2);
% % Tq = interp1(sorted_A2(I2),sorted_T2(I2),Aq);
% % plot(smooth(taxis(round(Tq)),5),smooth(phi(round(Aq)),5),'color','c','linewidth',1.5);hold on;

% xlim([27.5,27.5]+[0,30]);

set(gca,'xtick',[]);
%
% fig_current = gcf; fig_current.Renderer = 'painters';
% print(fig_current,['plot curvs/deltacurv_mincurvannotated_' '_cetD09_highthres'],'-dpdf');


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

function imArr = readfromtiff(filename,t_total,size_x,size_y)

imArr = zeros(size_x,size_y,t_total);

% mainfolder = cd(tiffdir);
for t = 1:1:t_total
    disp(t);
    I = imread([filename,'.tif'],t);
    imArr(:,:,t) = double(I);

end
% cd(mainfolder);

end
