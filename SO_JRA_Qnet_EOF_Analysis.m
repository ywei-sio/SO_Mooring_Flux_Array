clear;clc;
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


%% define analysis domain
LatBry=-20;
iy=find(lat<=LatBry);
lat0=lat(iy);
[XX0,YY0]=meshgrid(lon0,lat0);

% var=Qnet0(:,iy,:);

% load Qnet_HighFreq_Hanning  
% var=Qnet_HighFreq(:,iy,:);

% load Qnet_LowFreq
load Qnet_LowFreq_Hanning_smooth
var=Qnet_LowFreq(:,iy,:);


%---remove the mean---------
var_mean=mean(var,3);
for it=1:nt
   var(:,:,it)=var(:,:,it)-var_mean; 
end
%----------------------------
% flag=0;
% var_std=std(var,flag,3);
% for it=1:nt
%     var(:,:,it)=var(:,:,it)./var_std;
% end


%-----weight the variable ----------
W(1:nx,1:ny)=0; % weight for global field
for j=1:ny
    W(:,j)=cos(lat(j)/90*pi/2);
end
for it=1:nt
    var(:,:,it)=squeeze(var(:,:,it)).*W(:,iy);
end




[lonlen,latlen,timelen]=size(var);
aa=reshape(squeeze(var(:,:,1)),1,latlen*lonlen);% aa is a temporal var.
IX1=~isnan(aa); 
n=sum(IX1);

P(1:timelen,1:n)=nan;
for it=1:timelen
    aa=reshape(squeeze(var(:,:,it)),1,latlen*lonlen);% aa is a temporal var.
%     IX1=~isnan(aa);aa(isnan(aa))=[];
%     P(it,1:length(aa))=aa;   
P(it,1:n)=aa(IX1); 
end
% [m,n]=size(P);


[COEFF,SCORE,latent,tsquare] = princomp(P,'econ');%
frac=latent/sum(latent)*100;

%---------------------figure---------------------
figure('Position',[10 10 1000 800]);
cv=-7:.5:7;
lon_title=-28;    lat_title=-2;
A(1:lonlen,1:latlen,1:3)=nan;
%--------------------first mode-------------------------
x1=COEFF(:,1);bb=[];bb(1:latlen*lonlen)=nan;bb(IX1)=x1./std(x1,0);xx1=reshape(bb,lonlen,latlen);
xx1=xx1./W(:,iy); % divide the weight
y=SCORE(:,1);yy=y*std(x1,0);
A(1:lonlen,1:latlen,1)=xx1;

subplot('Position',[0.1 0.7 0.4 0.23]);
m_proj('stereographic','lat',-90,'long',0,'radius',70,'rec','off');hold on;
[~,h]=m_contourf(XX0,YY0,xx1',cv);set(h,'linestyle','none');colorbar;caxis([-5 5]);

m_grid('xtick',12,'XAxisLocation','top','tickdir','out','ytick',-80:20:-20,'linest','-','color','k');
m_coast('patch',[.7 .7 .7],'edgecolor','none');
set(findobj('tag','m_grid_color'),'facecolor','none');

 m_text(lon_title,lat_title,['1st EOF (Var:',num2str(frac(1),'%10.2f'),'%)'], 'fontsize',12,'fontweight','bold');
% m_text(118,40,'(a)','fontsize',12,'color','k');

subplot('Position',[0.55 0.72 0.4 0.2]);
plot((0:timelen-1)/12+1979+1/24,yy,'color',[0  0.4470  0.7410]);hold on
xlim([1979 2017]);ylim([-8 8]);set(gca,'ytick',-8:4:8)

hold on;plot([1979 2017],[0 0],'k:');
text(1980,7.0,'(d)','fontsize',12);

ylabel('(Wm^{-2})','fontsize',12)
xlabel('Time (year)','fontsize',12)


%--------------------second mode-------------------------
x1=COEFF(:,2);bb=[];bb(1:latlen*lonlen)=nan;bb(IX1)=x1./std(x1,0);xx1=reshape(bb,lonlen,latlen);
xx1=xx1./W(:,iy); % divide the weight

y=SCORE(:,2);yy=y*std(x1,0);
A(1:lonlen,1:latlen,2)=xx1;

subplot('Position',[0.1 0.4 0.4 0.23]);
m_proj('stereographic','lat',-90,'long',0,'radius',70,'rec','off');hold on;
[~,h]=m_contourf(XX0,YY0,xx1',cv);set(h,'linestyle','none');colorbar;caxis([-5 5]);

m_grid('xtick',12,'XAxisLocation','top','tickdir','out','ytick',-80:20:-20,'linest','-','color','k');
m_coast('patch',[.7 .7 .7],'edgecolor','none');
set(findobj('tag','m_grid_color'),'facecolor','none');
 m_text(lon_title,lat_title,['2nd EOF (Var:',num2str(frac(2),'%10.2f'),'%)'], 'fontsize',12,'fontweight','bold');

subplot('Position',[0.55 0.42 0.4 0.2]);
plot((0:timelen-1)/12+1979+1/24,yy,'color',[0  0.4470  0.7410]);hold on
xlim([1979 2017]);ylim([-8 8]);set(gca,'ytick',-8:4:8)

hold on;plot([1979 2017],[0 0],'k:');
text(1980,7.0,'(e)','fontsize',12);

ylabel('(Wm^{-2})','fontsize',12)
xlabel('Time (year)','fontsize',12)


%--------------------third mode-------------------------
x1=COEFF(:,3);bb=[];bb(1:latlen*lonlen)=nan;bb(IX1)=x1./std(x1,0);xx1=reshape(bb,lonlen,latlen);
xx1=xx1./W(:,iy); % divide the weight

y=SCORE(:,3);yy=y*std(x1,0);
A(1:lonlen,1:latlen,3)=xx1;


subplot('Position',[0.1 0.1 0.4 0.23]);
m_proj('stereographic','lat',-90,'long',0,'radius',70,'rec','off');hold on;
[c,h]=m_contourf(XX0,YY0,xx1',cv);set(h,'linestyle','none');colorbar;caxis([-5 5]);

m_grid('xtick',12,'XAxisLocation','top','tickdir','out','ytick',-80:20:-20,'linest','-','color','k');
m_coast('patch',[.7 .7 .7],'edgecolor','none');
set(findobj('tag','m_grid_color'),'facecolor','none');
m_text(lon_title,lat_title,['3rd EOF (Var:',num2str(frac(3),'%10.2f'),'%)'], 'fontsize',12,'fontweight','bold');


subplot('Position',[0.55 0.12 0.4 0.20]);
plot((0:timelen-1)/12+1979+1/24,yy,'color',[0  0.4470  0.7410]);hold on
xlim([1979 2017]);ylim([-8 8])
hold on;plot([1979 2017],[0 0],'k:');
text(1980,7.0,'(f)','fontsize',12);
set(gca,'ytick',-8:4:8)
ylabel('(Wm^{-2})','fontsize',12)
xlabel('Time (year)','fontsize',12)

