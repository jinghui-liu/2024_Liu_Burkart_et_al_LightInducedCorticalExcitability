global umperpix spf size_x size_y

umperpix = 303.64/512;
spf = 1;
size_x = 512;size_y = 512;

ROI = [51,100;255,285];

t_total = 1000;

%% Cry2-Null

% imArr = zeros(size_x,size_y,t_total);
% 
% rad=5; tfac=1; smfac=11; sFac=0.99; w=5; %%
% bws = cell(1,t_total);
% 
% for t = 1:1:t_total
%     
%     disp(t);
%     
%     I = imread(['Cry2_NULL-CIBN 31 4xdilute_M11 processed/pos2_actexp1_tseries1_lightat100_speed1_norepeat_spf1s_t1000-ch1.tif'],t);
%     
%     I = double(I);
%     
%     imArr(:,:,t) = I;
% 
%     rc1 = getboundary(I,rad,tfac,smfac,sFac);
%     
%     bws{t} = rc1;
%     
% %     figure;
% %     imagesc(I);daspect([1,1,1]);hold on;
% %     plot(rc1(:,2),rc1(:,1),'r-','linewidth',2);hold on;
% %     
%    
% end
% 
% 
% save('Cry2_NULL-CIBN 31 4xdilute_M11 processed/imArr_pos2_actexp1_tseries1_lightat100_speed1_norepeat_spf1s_t1000-ch1.mat',...
%     'imArr');
% 
% save('Cry2_NULL-CIBN 31 4xdilute_M11 processed/bws_pos2_actexp1_tseries1_lightat100_speed1_norepeat_spf1s_t1000-ch1.mat',...
%     'bws','rad','tfac','smfac','sFac');

%%

load('Cry2_NULL-CIBN 31 4xdilute_M11 processed/imArr_pos2_actexp1_tseries1_lightat100_speed1_norepeat_spf1s_t1000-ch1.mat',...
    'imArr');

load('Cry2_NULL-CIBN 31 4xdilute_M11 processed/bws_pos2_actexp1_tseries1_lightat100_speed1_norepeat_spf1s_t1000-ch1.mat',...
    'bws');


%%

t = 1;

clear bws_angle

% figure('position',[100,100,1000,1000]);
% imagesc(imArr(:,:,t));daspect([1,1,1]);hold on;
% 
% plot(bws{t}(:,2),bws{t}(:,1),'color','w','linewidth',1.5,'linestyle','--');hold on;
% 
c = mean(ROI(1,:));
r = mean(ROI(2,:));
% 
% % [c,r] =ginput(1);
% 
% scatter(c,r,100,'r');hold on;
% 
% D = pdist2(bws{t},[r,c]);
% minloc = find(D==min(D));
% 
% scatter(bws{t}(minloc,2),bws{t}(minloc,1),100,'r','filled');hold on;

pxlen = 300; % in pixel
subkymo = zeros(pxlen+1,t_total);
subkymo_oppo = zeros(pxlen+1,t_total);

subw = 10; % 10, 1, 8

phis_interpcor = (-38:1:38).*pi/180;
subkymo_interp = zeros(length(phis_interpcor),t_total);
subkymo_interp_norm = zeros(length(phis_interpcor),t_total);

for t = 1:1:t_total
    
    imArr_gauss = imgaussfilt(imArr(:,:,t),1);
    
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
    bws_angle{t} = (arclens./perimeter)*2*pi;
    
    
    
    
    
    D_oppo = pdist2(bws{t},2.*mean(bws{t})-[r,c]);
    minloc_oppo = find(D_oppo==min(D_oppo));minloc_oppo = minloc_oppo(1);
    
    if (minloc_oppo-pxlen/2)<1
        bws_kymo_oppo = [bws{t}(length(bws{t})+minloc_oppo-pxlen/2:end,:);bws{t}(1:minloc_oppo+pxlen/2,:)];
    else
        if (minloc_oppo+pxlen/2)>length(bws{t})
            bws_kymo_oppo = [bws{t}(minloc_oppo-pxlen/2:end,:);bws{t}(1:minloc_oppo+pxlen/2-length(bws{t}),:)];
        else
            bws_kymo_oppo = bws{t}(minloc_oppo-pxlen/2:minloc_oppo+pxlen/2,:);
        end
    end
    
    
    for i = 1:1:length(bws_kymo)

        subI = imArr_gauss(bws_kymo(i,1)-subw:bws_kymo(i,1)+subw,bws_kymo(i,2)-subw:bws_kymo(i,2)+subw);
        subI_sorted = sort(subI(:),'descend');
        
        subkymo(i,t) = mean(subI_sorted(1:round(length(subI_sorted)/8))); %/((2*subw+1)*(2*subw+1));
        
        subI_oppo = imArr_gauss(bws_kymo_oppo(i,1)-subw:bws_kymo_oppo(i,1)+subw,bws_kymo_oppo(i,2)-subw:bws_kymo_oppo(i,2)+subw);
        subI_sorted_oppo = sort(subI_oppo(:),'descend');
        
        subkymo_oppo(i,t) = mean(subI_sorted_oppo(1:round(length(subI_sorted_oppo)/8))); %/((2*subw+1)*(2*subw+1));
        
%         subkymo(i,t) = sum(sum(imArr(bws_kymo(i,1)-subw:bws_kymo(i,1)+subw,...
%             bws_kymo(i,2)-subw:bws_kymo(i,2)+subw,t)))./((2*subw+1)*(2*subw+1));
    end
    
    
    
    [uniq_phit,I_phi] = unique(bws_angle{t});
    
    subkymo_interp_t = interp1(uniq_phit,subkymo(I_phi,t),phis_interpcor);
    
    subkymo_interp(:,t) = subkymo_interp_t;
    
    subkymo_interp_norm(:,t) = (subkymo_interp(:,t)-ones(size(subkymo_interp(:,t))).*mean(subkymo_oppo(:,t)));



end

% figure;
% imagesc(imArr(:,:,end));hold on;daspect([1,1,1]);
% plot(bws_kymo(:,2),bws_kymo(:,1),'w-','linewidth',2);hold on;
% 
% 

%%

% smrad = 5;
% 
% mymap = colorcet('L06', 'N', 5);
% 
% figure;
% t = 100;
% h1 = plot(phis_interpcor,smooth(subkymo_interp_norm(:,t)./max(subkymo_interp_norm(:)),smrad),'linewidth',1.5,'color',mymap(1,:));hold on;
% t = 201;
% h2 = plot(phis_interpcor,smooth(subkymo_interp_norm(:,t)./max(subkymo_interp_norm(:)),smrad),'linewidth',1.5,'color',mymap(2,:));hold on;
% t = 401;
% h3 = plot(phis_interpcor,smooth(subkymo_interp_norm(:,t)./max(subkymo_interp_norm(:)),smrad),'linewidth',1.5,'color',mymap(3,:));hold on;
% t = 1000;
% h4 = plot(phis_interpcor,smooth(subkymo_interp_norm(:,t)./max(subkymo_interp_norm(:)),smrad),'linewidth',1.5,'color',mymap(4,:));hold on;
% 
% plot([-pi/4 pi/4],[0,0],'k--','linewidth',1,'handlevisibility','off');
% 
% xlim([-pi/4 pi/4]);ylim([-0.2 1]);
% set(gca,'xtick',[-pi/4 0 pi/4],'xticklabel',{'-pi/4','0','pi/4'});
% pbaspect([1,1,1]);
% 
% %
% 
% D = pdist2(bws{100},[r,c]);
% minloc = find(D==min(D));minloc = minloc(1);
% D = pdist2(bws{100},[285,51]);
% minloc1 = find(D==min(D));minloc1 = minloc1(1);
% D = pdist2(bws{100},[255,51]);
% minloc2 = find(D==min(D));minloc2 = minloc2(1);
% 
% bws_illuminated = [bws{100}(minloc1:minloc2,:)];
% 
% % figure;
% % imagesc(imArr(:,:,100));daspect([1,1,1]);hold on;
% % plot(bws_illuminated(:,2),bws_illuminated(:,1),'w-','linewidth',2);hold on;
% 
% perimeter = 0;
% for i = 1:1:length(bws{100})-1
%     perimeter = perimeter+pdist2(bws{100}(i,:),bws{100}(i+1,:));
% end
% perimeter = perimeter+pdist2(bws{100}(end,:),bws{100}(1,:));
%    
% clear arclens
% arclens(1) = 0;
% for i = 1:1:length(bws_illuminated)-1
%     arclens(i+1) = arclens(i)+pdist2(bws_illuminated(i,:),bws_illuminated(i+1,:));
% end
% arclens = arclens-arclens((minloc-minloc1+1));
% bws_angle = (arclens./perimeter)*2*pi;
% 
% plot([bws_angle(1),bws_angle(1)],[-0.2,1],'c--','linewidth',1,'handlevisibility','off');hold on;
% plot([bws_angle(end),bws_angle(end)],[-0.2,1],'c--','linewidth',1,'handlevisibility','off');hold on;
% 
% legend({'0s','100s','300s','900s'},'fontsize',14);
% 
% 
% 
% % fig_current = gcf; fig_current.Renderer = 'painters';
% % print(fig_current,['20221214_2b_OptoNull_expt'],'-dpdf');

%% Extended data figure, remove normalization

mymap = colorcet('L06', 'N', 5);

figure;
t = 100;
h1 = plot(phis_interpcor,subkymo_interp(:,t),'linewidth',1.5,'color',mymap(1,:));hold on;
t = 201;
h2 = plot(phis_interpcor,subkymo_interp(:,t),'linewidth',1.5,'color',mymap(2,:));hold on;
t = 401;
h3 = plot(phis_interpcor,subkymo_interp(:,t),'linewidth',1.5,'color',mymap(3,:));hold on;
t = 1000;
h4 = plot(phis_interpcor,subkymo_interp(:,t),'linewidth',1.5,'color',mymap(4,:));hold on;

% plot([-pi/4 pi/4],[0,0],'k--','linewidth',1,'handlevisibility','off');

xlim([-pi/4 pi/4]);% ylim([-0.2 1]);
set(gca,'xtick',[-pi/4 0 pi/4],'xticklabel',{'-pi/4','0','pi/4'});
pbaspect([1,1,1]);



D = pdist2(bws{100},[r,c]);
minloc = find(D==min(D));minloc = minloc(1);
D = pdist2(bws{100},[285,51]);
minloc1 = find(D==min(D));minloc1 = minloc1(1);
D = pdist2(bws{100},[255,51]);
minloc2 = find(D==min(D));minloc2 = minloc2(1);

bws_illuminated = [bws{100}(minloc1:minloc2,:)];

% figure;
% imagesc(imArr(:,:,100));daspect([1,1,1]);hold on;
% plot(bws_illuminated(:,2),bws_illuminated(:,1),'w-','linewidth',2);hold on;

perimeter = 0;
for i = 1:1:length(bws{100})-1
    perimeter = perimeter+pdist2(bws{100}(i,:),bws{100}(i+1,:));
end
perimeter = perimeter+pdist2(bws{100}(end,:),bws{100}(1,:));
   
clear arclens
arclens(1) = 0;
for i = 1:1:length(bws_illuminated)-1
    arclens(i+1) = arclens(i)+pdist2(bws_illuminated(i,:),bws_illuminated(i+1,:));
end
arclens = arclens-arclens((minloc-minloc1+1));
bws_angle = (arclens./perimeter)*2*pi;

plot([bws_angle(1),bws_angle(1)],[10,20],'c--','linewidth',1,'handlevisibility','off');hold on;
plot([bws_angle(end),bws_angle(end)],[10,20],'c--','linewidth',1,'handlevisibility','off');hold on;

legend({'0s','100s','300s','900s'},'fontsize',14);

% fig_current = gcf; fig_current.Renderer = 'painters';
% print(fig_current,['20230106_S4e_OptoNull_expt_unnormalized'],'-dpdf');
