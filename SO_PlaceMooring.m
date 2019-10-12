clear;clc;
close all;
%-------------define some variables-------------------------------------
LatBry=-20; % set this value smaller than -20
Do_EEZ_Mask=1; % mask EEZ (and shallow water) or not

De_OOI=0;
De_SOFS=0;% remove OOI and SOFS site correlated component or not

Freqcomp='Low';% 'High' freq. or 'Low' freq.
Choicerule='Totvar'; % site mooring to maximize local SD or total variance. Choicerule='SD' or 'TotVar'

%------ mooring selection experiment based on JRA reanalysis------
load JRA_Qnet_197901_201612.mat;
nt=456;
nx=length(lon);
ny=length(lat);
nyear=2016-1979+1;

missval=-1.e34;
landmask(landmask==0)=nan;
Qnet(nx,ny,nt)=0;
for it=1:nt
    Qnet(:,:,it)=squeeze(Qnet(:,:,it)).*landmask;
end
Qnet0=Qnet;Qnet0(1:nx/2,:,:)=Qnet(nx/2+1:nx,:,:);Qnet0(nx/2+1:nx,:,:)=Qnet(1:nx/2,:,:);
lon0=lon;lon0(1:nx/2)=lon(nx/2+1:nx)-360;lon0(nx/2+1:nx)=lon(1:nx/2);

lon0(end)=lon0(end)+0.5; % make it close to 180 to avoid space in figure
landmask0=landmask;landmask0(1:nx/2,:)=landmask(nx/2+1:nx,:);landmask0(nx/2+1:nx,:)=landmask(1:nx/2,:);

%% define three basins
%-------not used so far-------------
lon1=148;lon2=-70.9;lon3=22;
ix_pac=find(lon0>lon1|lon0<lon2);
ix_alt=find(lon0>lon2&lon0<lon3);
ix_ind=find(lon0>lon3&lon0<lon1);
%% define analysis domain
% LatBry=-20;
iy=find(lat<=LatBry);
lat0=lat(iy);
[XX0,YY0]=meshgrid(lon0,lat0);

%% mask EEZ and shallow water (less than 2000 m)
load  JRA_EEZ_Topo_Mask %EEZ_Topo_Mask, south of 20S
EEZ_Topo_Mask0=EEZ_Topo_Mask(:,iy);
%% define the possible domain to deploy moorings
iy_band=find(lat>=-35|lat<=-65);
JRA_mooring_band(1:nx,1:ny)=0;JRA_mooring_band(:,iy_band)=1;
Mooring_band=JRA_mooring_band(:,iy);


%% load data
if strcmp(Freqcomp,'High')==1
    load Qnet_HighFreq_Hanning;
    var_global=Qnet_HighFreq;
    %variables for plot
    cv_sd=0:.5:30; cmin_sd=0; cmax_sd=40;
    cv_var=0:1:50; cmin_var=0;cmax_var=30;
    cv_corr=-1:.1:1;

elseif strcmp(Freqcomp,'Low')==1
    load Qnet_LowFreq_Hanning_smooth;
    var_global=Qnet_LowFreq;
    
    %variables for plot
    cv_sd=0:.5:30; cmin_sd=0; cmax_sd=10;
    cv_var=0:1:50; cmin_var=0;cmax_var=10;
    cv_corr=-1:.1:1;

end
var=var_global(:,iy,:); % variable in the Southern Ocean


%% ---------the position of OOI and SOFS moorings-------------
ooi_lon=-89.28;ooi_lat=-54.47;
sofs_lat=-47; sofs_lon=142;

%----Interpolate the data at ooi and sofs sites and get their correlated components ---------------
var_ooi(1:nt)=missval;
var_sofs(1:nt)=missval;

for i=1:nt
    var_ooi(i)=interp2(lon0,lat,squeeze(var_global(:,:,i))',ooi_lon,ooi_lat);
    var_sofs(i)=interp2(lon0,lat,squeeze(var_global(:,:,i))',sofs_lon,sofs_lat);
end

reg_ooi(1:nx,1:ny,1:nt)=nan;
reg_sofs(1:nx,1:ny,1:nt)=nan;
reg_ooi_sofs(1:nx,1:ny,1:nt)=nan;

% warning off;
for i=1:nx
    for j=1:ny
        if ~isnan(landmask0(i,j))
            [b,~] = regress(squeeze(var_global(i,j,:)),[var_ooi',ones(nt,1)]);
            reg_ooi(i,j,:)=b(1)*var_ooi+b(2)*ones(nt,1)';
            
            [b,~] = regress(squeeze(var_global(i,j,:)),[var_sofs',ones(nt,1)]);
            reg_sofs(i,j,:)=b(1)*var_sofs+b(2)*ones(nt,1)';
            
            [b,~] = regress(squeeze(var_global(i,j,:)),[var_ooi',var_sofs',ones(nt,1)]);
            reg_ooi_sofs(i,j,:)=b(1)*var_ooi+b(2)*var_sofs+b(3)*ones(nt,1)';
        end
    end
end



flag=0;% std flag
if De_OOI==1&&De_SOFS==0
    var=var-reg_ooi(:,iy,:);
end
if De_SOFS==1&&De_OOI==0
    var=var-reg_sofs(:,iy,:);
end
if De_SOFS==1&&De_OOI==1
    var=var-reg_ooi_sofs(:,iy,:);
end


%% --------remove the mean------------
var_mean=mean(var,3);
for it=1:nt
    var(:,:,it)=var(:,:,it)-var_mean;
end
%----------------------------------

%
% var_std=std(var,flag,3);
% for it=1:nt
%     var(:,:,it)=var(:,:,it)./var_std;
% end
%

%%
[lonlen,latlen,timelen]=size(var);
IX1=isnan(var);% 3-dimension

aa=reshape(squeeze(var(:,:,1)),1,latlen*lonlen);% aa is a temporal var.
IX2=~isnan(aa);
n=sum(IX2);

P(1:timelen,1:n)=nan;
for it=1:timelen
    aa=reshape(squeeze(var(:,:,it)),1,latlen*lonlen);% aa is a temporal var.
    P(it,1:n)=aa(IX2);
end
[m,n]=size(P); %m=timelen, n denotes grid number in the ocean




var_reg(1:lonlen,1:latlen,1:nt)=nan;% mooring site-correlated component that neeeds to be subtracted
var_corr(1:lonlen,1:latlen)=nan; % correlation map

% Weight(1:lonlen,1:latlen)=0;
Var_total(1:lonlen,1:latlen)=nan;% total variance can be explained by each grid point
Var_ave(1:lonlen,1:latlen)=nan;% total variance can be explained by each grid point

%%
dplot(1:nx,1:ny)=nan;

Lx=5;Ly=5;
dx=40;dy=40;

it_num=8; % iteration number or mooring number
cord_lon_lat(1:it_num,1:2)=0;
var_sum(1:it_num)=0;
%% keep distance from selected moorings,
Mooring_Cover(1:lonlen,1:latlen)=0; % the area where physical processes can be resolved with mooring sites, 
Mooring_Distance(1:lonlen,1:latlen)=0; % the distance of the grid points to a single mooring site



%load cmap_pos_neg
%load camp_positive

% draw different variables in 3 figures with different colormaps  
hf1=figure('position',[100 50 900 950]);
hf2=figure('position',[100 50 900 950]);
hf3=figure('position',[100 50 900 950]);

A(1:lonlen,1:latlen,1:3,1:it_num)=nan; % used to store variables generated in iteration process
for iter=1:it_num
    %     plot1: std of var，and maximum value of std
    %     plot2: domain averaged variance explained by mooring at each grid
    %     plot3: correlation map
    for it=1:timelen
        aa=reshape(squeeze(var(:,:,it)),1,latlen*lonlen);% aa is a temporal var.
        P(it,1:n)=aa(IX2);
    end
    
    
    for id=1:lonlen
        for jd=1:latlen
            xtmp=squeeze(var(id,jd,:));
            if ~isnan(xtmp(1))
%                 sum=0;
%                 for k=1:n
%                     sum=sum+(dot(P(:,k),xtmp))^2/(dot(xtmp,xtmp));
%                 end
%                 sum/m/n
                
                xtmp=xtmp/norm(xtmp);
                reg_coef=(P)'*xtmp;  % reg_coef of SO Qnet field to this site
                Var_total(id,jd)=norm(reg_coef).^2/m; % m is length of time
                Var_ave(id,jd)=Var_total(id,jd)/n;% n is grid point number; make it stands for grid point-averaged variance.
                
            end
        end
    end
    
    dplot=std(var,flag,3); % calculate standard deviation of Qnet

%------------update the mask-------------------------- 
%1. land or sea % isnan(dplot)
%2. EEZ/shallow water %  EEZ_Topo_Mask
%3. Domain to deploy mooring % Mooring_band
%4. Keep a distance from existing mooring sites % Mooring_Select   
 
if Do_EEZ_Mask==1
    IX20=isnan(dplot)|EEZ_Topo_Mask0|Mooring_band|Mooring_Cover;
else
    % IX20=isnan(dplot);
    IX20=isnan(dplot)|Mooring_band|Mooring_Cover;
end

if strcmp(Choicerule,'SD')==1
    dplotS=dplot; % set it to SD
    % dplotS=smooth2D_per(dplot,Lx,Ly,lat0);
elseif strcmp(Choicerule,'Totvar')==1
    dplotS=Var_ave; % set it to overall variance explained
else
    dplotS=Var_ave; % set it to overall variance explained
end

dplotS(IX20)=missval;% change nan to missval so that we can locate the position of max


[var_max1,ix_max1]=max(dplotS);
[var_max2,ix_max2]=max(var_max1);
% locate maximum value
dy0=lat0(ix_max2);
dx0=lon0(ix_max1(ix_max2));
var_max=squeeze(var(ix_max1(ix_max2),ix_max2,:));


%record the coordinate of mooring sites
cord_lon_lat(iter,1)=dx0; cord_lon_lat(iter,2)=dy0;
    
% calcualte the correlated component and corrleation coefficient to
% variable at this mooring site
for i=1:lonlen
    for j=1:latlen
        if ~isnan(dplot(i,j))
            [b,bint]=regress(squeeze(var(i,j,:)),var_max);
            var_reg(i,j,:)=b*var_max;
            
            r=corrcoef(squeeze(var(i,j,:)),var_max);
            var_corr(i,j)=r(1,2);
        end
    end
end


%     %--------add on a weight-----------
%     i=ix_max1(ix_max2);
%     j=ix_max2;
%     Weight(1:lonlen,1:latlen)=0;
%     ix1_w=mod(i-dx-1:i-1,lonlen)+1;ix2_w=mod(i:i+dx-1,lonlen)+1;
%     iy1_w=max(1,j-dy);iy2_w=min(latlen,j+dy);
% %     Wx=exp(-abs(-dx:1:dx)/dx);Wy=exp(-abs((iy1_w:iy2_w)-j)/dy);
% %     Weight([ix1_w,ix2_w],iy1_w:iy2_w)=Wx'*Wy;
%     Weight([ix1_w,ix2_w],iy1_w:iy2_w)=1;
%     for it=1:nt
%          var_reg(:,:,it)=squeeze(var_reg(:,:,it)).*Weight;
%     end


%%
figure(hf1)
subplot(it_num,3,(iter-1)*3+1);

m_proj('stereographic','lat',-90,'long',0,'radius',LatBry+90,'rec','off');hold on;
[~,h]=m_contourf(XX0,YY0,dplot',cv_sd);set(h,'linestyle','none');colorbar;caxis([cmin_sd cmax_sd]);
colormap cool;

m_grid('xtick',12,'XAxisLocation','top','tickdir','out','ytick',-80:20:-20,'yticklabel',[],'linest','-','color','k');
m_coast('patch',[.7 .7 .7],'edgecolor','none');
set(findobj('tag','m_grid_color'),'facecolor','none');

m_plot(dx0,dy0,'k-o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5)
    
%%
figure(hf2)
subplot(it_num,3,(iter-1)*3+2);
m_proj('stereographic','lat',-90,'long',0,'radius',LatBry+90,'rec','off');hold on;
[~,h]=m_contourf(lon0,lat0,Var_ave',cv_var);set(h,'linestyle','none');caxis([cmin_var cmax_var]);colorbar

m_grid('xtick',12,'XAxisLocation','top','tickdir','out','ytick',-80:20:-20,'yticklabel',[],'linest','-','color','k');
m_coast('patch',[.7 .7 .7],'edgecolor','none');
set(findobj('tag','m_grid_color'),'facecolor','none');

colormap cool;
m_plot(dx0,dy0,'k-o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5)

%%
figure(hf3)
subplot(it_num,3,(iter-1)*3+3);
m_proj('stereographic','lat',-90,'long',0,'radius',LatBry+90,'rec','off');hold on;
[c,h]=m_contourf(lon0,lat0,var_corr',cv_corr);set(h,'linestyle','none');caxis([-1 1]);colorbar

m_grid('xtick',12,'XAxisLocation','top','tickdir','out','ytick',-80:20:-20,'yticklabel',[],'linest','-','color','k');
m_coast('patch',[.7 .7 .7],'edgecolor','none');
set(findobj('tag','m_grid_color'),'facecolor','none');
colormap jet;
m_plot(dx0,dy0,'k-o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5)
%% update the variable for next iteration   
var=var-var_reg;
var=var+0.0000000001*rand(lonlen,latlen,nt);% added on a small random variable to avoid zero denominator 
var(IX1)=nan;

%mask the area around the selected moorings

%option 1: keep a 15 longitude degree from mooring  
i=ix_max1(ix_max2);
dx=27; % keep a 15 longitude degree from mooring, which approximates 27 grid points
ix1_w=mod(i-dx-1:i-1,lonlen)+1;ix2_w=mod(i:i+dx-1,lonlen)+1;
Mooring_Cover([ix1_w,ix2_w],:)=1;

% %option 2: keep a distance from mooring site, e.g., 500km
% for i=1:lonlen
%     for j=1:latlen
%         [dist,phaseangle] = sw_dist([lat0(j),dy0],[lon0(i),dx0], 'km');
%         Mooring_Distance(i,j)=dist;
%     end
% end
% Mooring_Cover(Mooring_Distance<500)=1; 
%%
A(:,:,1,iter)= dplot;
A(:,:,2,iter)= Var_ave;
A(:,:,3,iter)= var_corr;
end

save cord_lon_lat_TotVar_20S_HighFreq  cord_lon_lat


% set(hf1,'position',[100 50 900 950]);
figure;
m_proj('stereographic','lat',-90,'long',0,'radius',70,'rec','off');hold on;
m_grid('xtick',12,'XAxisLocation','top','tickdir','out','ytick',-80:20:-20,'yticklabel',[],'linest','-','color','k');
m_coast('patch',[.7 .7 .7],'edgecolor','none');

% h1=m_plot(cord_lon_lat(1,1),cord_lon_lat(1,2),'o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',10);
% h2=m_plot(cord_lon_lat(2,1),cord_lon_lat(2,2),'p','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',10);
h1=m_plot(cord_lon_lat(1,1),cord_lon_lat(1,2),'o','MarkerEdgeColor','none','MarkerFaceColor','g','MarkerSize',10);
h2=m_plot(cord_lon_lat(2,1),cord_lon_lat(2,2),'p','MarkerEdgeColor','none','MarkerFaceColor','g','MarkerSize',10);
h3=m_plot(cord_lon_lat(3,1),cord_lon_lat(3,2),'o','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',10);
h4=m_plot(cord_lon_lat(4,1),cord_lon_lat(4,2),'p','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',10);
h5=m_plot(cord_lon_lat(5,1),cord_lon_lat(5,2),'p','MarkerEdgeColor','none','MarkerFaceColor',[160 32 240]/250,'MarkerSize',10);

h=legend([h1,h2,h3,h4,h5],'1','2','3','4','5');
% set(h, 'Box', 'off','Orientation','horizontal','Location','Northeast');




