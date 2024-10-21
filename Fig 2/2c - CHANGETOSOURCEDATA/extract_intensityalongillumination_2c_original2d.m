
global umperpix spf

umperpix = 303.64/512;
spf = 1;
size_x = 512;size_y = 512;

ROI = [410,459;243 273];

t_total = 1000;

filename = 'cry2ect2_M11';

%% Cry2-Ect2

% imArr = zeros(size_x,size_y,t_total);
% 
% rad=5; tfac=1; smfac=11; sFac=0.99; w=5; %%
% bws = cell(1,t_total);

% for t = 1:1:t_total
%     disp(t);
%     I = imread(['Ect2-CIBN 31 4xdilute_M11 processed/pos6_actexp1_tseries2_lightat1_speed1_norepeat_spf1s_t1000-ch1.tif'],t);
%     imArr(:,:,t) = double(I);
%     
%     rc1 = getboundary(I,rad,tfac,smfac,sFac);
    
%     bws{t} = rc1;
    
%     figure;
%     imagesc(I);daspect([1,1,1]);hold on;
%     plot(rc1(:,2),rc1(:,1),'r-','linewidth',2);hold on;
%   

% end

% save('Ect2-CIBN 31 4xdilute_M11 processed/imArr_pos6_actexp1_tseries2_lightat1_speed1_norepeat_spf1s_t1000-ch1.mat',...
%     'imArr');

% save('Ect2-CIBN 31 4xdilute_M11 processed/bws_pos6_actexp1_tseries2_lightat1_speed1_norepeat_spf1s_t1000-ch1.mat',...
%     'bws','rad','tfac','smfac','sFac');

%%

load('Ect2-CIBN 31 4xdilute_M11 processed/imArr_pos6_actexp1_tseries2_lightat1_speed1_norepeat_spf1s_t1000-ch1.mat',...
    'imArr');

load('Ect2-CIBN 31 4xdilute_M11 processed/bws_pos6_actexp1_tseries2_lightat1_speed1_norepeat_spf1s_t1000-ch1.mat',...
    'bws','rad','tfac','smfac','sFac');

%%

t = 1;

clear bws_angle

figure('position',[100,100,1000,1000]);
imagesc(imArr(:,:,t));daspect([1,1,1]);hold on;

plot(bws{t}(:,2),bws{t}(:,1),'color','w','linewidth',1.5,'linestyle','--');hold on;
% 
c = mean(ROI(1,:));
r = mean(ROI(2,:));

plot([ROI(1,1),ROI(1,2)],[ROI(2,1),ROI(2,1)],'c-','linewidth',2);hold on;
plot([ROI(1,1),ROI(1,2)],[ROI(2,2),ROI(2,2)],'c-','linewidth',2);hold on;
% [c,r] =ginput(1);

scatter(c,r,100,'r');hold on;

D = pdist2(bws{t},[r,c]);
minloc = find(D==min(D));

scatter(bws{t}(minloc,2),bws{t}(minloc,1),100,'r','filled');hold on;

%%

pxlen = 300; % in pixel
subkymo = zeros(pxlen+1,t_total);

subw = 10; % 10, 1, 8

phis_interpcor = (-38:1:38).*pi/180;
mindist_interp = zeros(length(phis_interpcor),t_total);
% subkymo_interp_norm = zeros(length(phis_interpcor),t_total);

for t = 1:1:t_total
    
    disp(t);
    
    D = pdist2(bws{t},[r,c]);
    minloc = find(D==min(D));minloc = minloc(1);
    
    if (minloc-pxlen/2)<1
        bws_kymo = [bws{t}(length(bws{t})+minloc-pxlen/2:end,:);bws{t}(1:minloc+pxlen/2,:)];
    else
        if (minloc+pxlen/2)>length(bws{t})
            bws_kymo = [bws{t}(minloc-pxlen/2:end,:);bws{t}(1:minloc+pxlen/2-length(bws{t}),:)];
        else
            bws_kymo = bws{t}(minloc-pxlen/2:minloc+pxlen/2,:);
        end
    end
    
    if t==1
        bws_kymo0 = bws_kymo;
        perimeter = 0;
        for i = 1:1:length(bws{t})-1
            perimeter = perimeter+pdist2(bws{t}(i,:),bws{t}(i+1,:));
        end
        perimeter = perimeter+pdist2(bws{t}(end,:),bws{t}(1,:));
        
        arclens(1) = 0;
        for i = 1:1:length(bws_kymo)-1
            arclens(i+1) = arclens(i)+pdist2(bws_kymo(i,:),bws_kymo(i+1,:));
        end
        arclens = arclens-arclens(pxlen/2+1);
        bws_angle0 = (arclens./perimeter)*2*pi;
    end
    
    for i = 1:1:length(bws_kymo0)
        D = pdist2(bws_kymo,bws_kymo0(i,:));
        minloc = find(D==min(D));minloc = minloc(1);
        
        mindist(i,t) = (pdist2(bws_kymo(minloc,:),bws_kymo0(i,:)))*umperpix;
        
    end
    
    [uniq_phit,I_phi] = unique(bws_angle0);
    
    mindist_interp_t = interp1(uniq_phit,mindist(I_phi,t),phis_interpcor);
    
    mindist_interp(:,t) = mindist_interp_t;
    
    
end

%%

% smrad = 5;
% 
% mymap = colorcet('L06', 'N', 5);
% 
% figure;
% t = 1;
% h1 = plot(phis_interpcor,smooth(mindist_interp(:,t)./max(mindist_interp(:)),smrad),'linewidth',1.5,'color',mymap(1,:));hold on;
% t = 101;
% h2 = plot(phis_interpcor,smooth(mindist_interp(:,t)./max(mindist_interp(:)),smrad),'linewidth',1.5,'color',mymap(2,:));hold on;
% t = 301;
% h3 = plot(phis_interpcor,smooth(mindist_interp(:,t)./max(mindist_interp(:)),smrad),'linewidth',1.5,'color',mymap(3,:));hold on;
% t = 901;
% h4 = plot(phis_interpcor,smooth(mindist_interp(:,t)./max(mindist_interp(:)),smrad),'linewidth',1.5,'color',mymap(4,:));hold on;
% 
% plot([-pi/4 pi/4],[0,0],'k--','linewidth',1,'handlevisibility','off');
% 
% xlim([-pi/4 pi/4]);ylim([-0.2 1]);
% set(gca,'xtick',[-pi/4 0 pi/4],'xticklabel',{'-pi/4','0','pi/4'});
% pbaspect([1,1,1]);
% 
% 
% D = pdist2(bws{1},[r,c]);
% minloc = find(D==min(D));minloc = minloc(1);
% D = pdist2(bws{1},[243,459]);
% minloc1 = find(D==min(D));minloc1 = minloc1(1);
% D = pdist2(bws{1},[273,459]);
% minloc2 = find(D==min(D));minloc2 = minloc2(1);
% 
% bws_illuminated = [bws{1}(minloc1:minloc2,:)];
% 
% % figure;
% % imagesc(imArr(:,:,1));daspect([1,1,1]);hold on;
% % plot(bws_illuminated(:,2),bws_illuminated(:,1),'w-','linewidth',2);hold on;
% 
% perimeter = 0;
% for i = 1:1:length(bws{1})-1
%     perimeter = perimeter+pdist2(bws{1}(i,:),bws{1}(i+1,:));
% end
% perimeter = perimeter+pdist2(bws{1}(end,:),bws{1}(1,:));
%    
% clear arclens
% arclens(1) = 0;
% for i = 1:1:length(bws_illuminated)-1
%     arclens(i+1) = arclens(i)+pdist2(bws_illuminated(i,:),bws_illuminated(i+1,:));
% end
% arclens = arclens-arclens((minloc-minloc1+1));
% bws_angle_illuminate = (arclens./perimeter)*2*pi;
% 
% plot([bws_angle_illuminate(1),bws_angle_illuminate(1)],[-0.2,1],'c--','linewidth',1,'handlevisibility','off');hold on;
% plot([bws_angle_illuminate(end),bws_angle_illuminate(end)],[-0.2,1],'c--','linewidth',1,'handlevisibility','off');hold on;
% 
% legend({'0s','100s','300s','900s'},'fontsize',14);
% 
% % fig_current = gcf; fig_current.Renderer = 'painters';
% % print(fig_current,['20221214_2d_OptoEct2_expt'],'-dpdf');

%% unnormalized for extended data

mymap = colorcet('L06', 'N', 5);

figure;
t = 1;
h1 = plot(phis_interpcor,mindist_interp(:,t),'linewidth',1.5,'color',mymap(1,:));hold on;
t = 101;
h2 = plot(phis_interpcor,mindist_interp(:,t),'linewidth',1.5,'color',mymap(2,:));hold on;
t = 301;
h3 = plot(phis_interpcor,mindist_interp(:,t),'linewidth',1.5,'color',mymap(3,:));hold on;
t = 901;
h4 = plot(phis_interpcor,mindist_interp(:,t),'linewidth',1.5,'color',mymap(4,:));hold on;

plot([-pi/4 pi/4],[0,0],'k--','linewidth',1,'handlevisibility','off');

xlim([-pi/4 pi/4]);% ylim([-0.2 1]);
set(gca,'xtick',[-pi/4 0 pi/4],'xticklabel',{'-pi/4','0','pi/4'});
pbaspect([1,1,1]);


D = pdist2(bws{1},[r,c]);
minloc = find(D==min(D));minloc = minloc(1);
D = pdist2(bws{1},[243,459]);
minloc1 = find(D==min(D));minloc1 = minloc1(1);
D = pdist2(bws{1},[273,459]);
minloc2 = find(D==min(D));minloc2 = minloc2(1);

bws_illuminated = [bws{1}(minloc1:minloc2,:)];

% figure;
% imagesc(imArr(:,:,1));daspect([1,1,1]);hold on;
% plot(bws_illuminated(:,2),bws_illuminated(:,1),'w-','linewidth',2);hold on;

perimeter = 0;
for i = 1:1:length(bws{1})-1
    perimeter = perimeter+pdist2(bws{1}(i,:),bws{1}(i+1,:));
end
perimeter = perimeter+pdist2(bws{1}(end,:),bws{1}(1,:));
   
clear arclens
arclens(1) = 0;
for i = 1:1:length(bws_illuminated)-1
    arclens(i+1) = arclens(i)+pdist2(bws_illuminated(i,:),bws_illuminated(i+1,:));
end
arclens = arclens-arclens((minloc-minloc1+1));
bws_angle_illuminate = (arclens./perimeter)*2*pi;

% plot([bws_angle_illuminate(1),bws_angle_illuminate(1)],[-0.2,1],'c--','linewidth',1,'handlevisibility','off');hold on;
% plot([bws_angle_illuminate(end),bws_angle_illuminate(end)],[-0.2,1],'c--','linewidth',1,'handlevisibility','off');hold on;

legend({'0s','100s','300s','900s'},'fontsize',14);

ylim([-1 8]);

% fig_current = gcf; fig_current.Renderer = 'painters';
% print(fig_current,['20230106_S4e_OptoEct2_expt_unnormalized'],'-dpdf');
