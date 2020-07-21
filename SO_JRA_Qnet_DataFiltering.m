clear;clc;
close all;
%----smooth the data--------- 
load JRA_Qnet_197901_201612.mat;

nx=length(lon);
ny=length(lat);
nyear=2016-1979+1;
nt=nyear*12;% 456


missval=-1.e34;
landmask(landmask==0)=nan;

% Qnet(nx,ny,nt)=0;
for it=1:nt
    Qnet(:,:,it)=squeeze(Qnet(:,:,it)).*landmask;
end
%----rearrange the data for plot-----------
Qnet0=Qnet;Qnet0(1:nx/2,:,:)=Qnet(nx/2+1:nx,:,:);Qnet0(nx/2+1:nx,:,:)=Qnet(1:nx/2,:,:);
lon0=lon;lon0(1:nx/2)=lon(nx/2+1:nx)-360;lon0(nx/2+1:nx)=lon(1:nx/2);
lon0(end)=lon0(end)+0.5; % make it close to 180 to avoid space in figure


%-----set the analysis domain-------------
LatBry=-20;
iy=find(lat<=LatBry);
lat0=lat(iy);
[XX0,YY0]=meshgrid(lon0,lat0);

%% ---------test filter based on OOI and SOFS moorings----------- 
% ---------positions of moorings and time of observation-------------
ooi_lon=-89.28;ooi_lat=-54.47;
sofs_lat=-47; sofs_lon=142;


t=(0:456-1)/12.0+1979+1/24;
ooi_date=(0:46-1)/12.0+2015+2/12+1/24;
sofs_date=(0:96-1)/12.0+2010+1/24;
%------------- the parameter for gft------------------
pvar = 1;dvar= .01;ifwd = 1;
kx1 = (2*pi)*1;%/365.2422;%wavenumbers to fit to: 1 yr   
kx2 = (2*pi)*2;%/365.2422;%wavenumbers to fit to: 6 mnt
kx3 = (2*pi)*3;%/365.2422;%wavenumbers to fit to: 4 mnt
kx4 = (2*pi)*4;%/365.2422;%wavenumbers to fit to: 3 mnt
% kx=[kx1 kx2 kx3 kx4]';
kx=[kx1 kx2]';
% kx=[kx1]';

%----------------interpolate to the mooring sites--------------------------
JRA_Qnet_ooi(1:456)=nan;JRA_Qnet_sofs(1:456)=nan;
nt=456;%
for i=1:nt
    JRA_Qnet_ooi(i)=interp2(lon0,lat,squeeze(Qnet0(:,:,i))',ooi_lon,ooi_lat);
    JRA_Qnet_sofs(i)=interp2(lon0,lat,squeeze(Qnet0(:,:,i))',sofs_lon,sofs_lat);
end
%----------------------remove annual signal--------------------------------
tmp1 = JRA_Qnet_ooi;
% tmp1=JRA_Qnet_sofs;
% tmp1a = detrend(tmp1 - mean(tmp1));
tmp1a = tmp1-mean(tmp1);

[p1,y1] = gft(t',tmp1a',kx,pvar,dvar,ifwd);
tmp1a = tmp1a - y1';% remove annual signals

%% --------------Low-pass filter OOI data and plot----------------
%-----Mirror extending the signal to suppress end-point problem----
L=length(tmp1a);
tmp1a0(1:L)=flip(tmp1a);tmp1a0(L+1:2*L)=tmp1a;tmp1a0(2*L+1:3*L)=flip(tmp1a); 

f_width=19;
tmp2a0=filter(hanning(f_width)/sum(hanning(f_width)),1,tmp1a0);
tmp2b0=tmp2a0(L+(f_width-1)/2+1:2*L+(f_width-1)/2);

% ------------------ figure---------------------------
figure('Position',[100,100,1200,600]);
subplot('Position',[.1 .6 .6 .35]);
plot(t,tmp1,'c-o','MarkerEdgeColor','c','MarkerFaceColor','c','MarkerSize',3);hold on
plot(t,y1,'r-o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',3);hold on
plot(t,tmp1a,'b-o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',3);hold on

ylabel('W m^{-2}','fontsize',12);
xlim([1978 2020]);
ylim([-250 250]);
text(1980,-200,'a','FontWeight','bold','fontsize',12);
set(gca,'TickDir','out','XMinorTick','on','YMinorTick','on','fontsize',12);
h=legend('Orginal','Harmonic','Residual');
set(h, 'Box', 'off','Orientation','horizontal','Location','Northeast','fontsize',10);

subplot('Position',[.10 .15 .6 .35]);
plot(t,tmp1a,'b-o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',3);hold on
plot(t,tmp2b0,'r','linewidth',2);
xlim([1978 2020]);
ylim([-105 105]);

ylabel('W m^{-2}','fontsize',12);
xlabel('Time (year)','fontsize',12)

text(1980,-80,'b','FontWeight','bold','fontsize',12);
set(gca,'TickDir','out','XMinorTick','on','YMinorTick','on','fontsize',12);%grid on

h=legend('Residual','LowFreq.');set(h, 'Box', 'off','Orientation','horizontal','Location','Northeast','fontsize',10);

%----- Burg power spectrum----------------------
subplot('Position',[.78 .15 .18 .35]);
box on

x1=detrend(tmp1a);x2=detrend(tmp2b0);

[Pxx,F,Pxxc] = pburg(x1,60,256,12,'ConfidenceLevel',0.90);

plot(F,10*log10(Pxxc(:,1)),'b:');hold on
plot(F,10*log10(Pxxc(:,2)),'b:');
h1=plot(F,10*log10(Pxx),'linewidth',1);


[Pxx,F,Pxxc] = pburg(x2,60,256,12,'ConfidenceLevel',0.90);

h2=plot(F,10*log10(Pxx),'r','linewidth',1);hold on
plot(F,10*log10(Pxxc(:,1)),'r:');
plot(F,10*log10(Pxxc(:,2)),'r:');

h=legend([h1,h2],'Residual','LowFreq');
set(h, 'Box', 'off','Location','Northeast','Orientation','horizontal','fontsize',10);

xlabel('Frequency/ cpy','fontsize',11)
ylabel('PSD/ 10 \times log_{10} (W m^{-2}/cpy)^{2}','fontsize',11)

set(gca,'ytick',-100:20:100,'ylim',[-100 50],'xlim',[0 10],'tickdir','out','fontsize',10);
set(gca,'XScale','log');

text(0.08,-80,'c','FontWeight','bold','fontsize',12);

%% ---------Processing Qnet data in the Southern Ocean------------------------
% ----Subtracting harmonic components to obtain HighFreq and then applying a hanning low-pass filter to HighFreq to get the LowFreq ------------------------------

[nx,ny,nt]=size(Qnet0);
Qnet_HighFreq(1:nx,1:ny,1:nt)=nan;
Qnet_LowFreq(1:nx,1:ny,1:nt)=nan;
Qnet_Harmonic(1:nx,1:ny,1:nt)=nan;

for i=1:nx
    for j=1:ny
        if ~isnan(Qnet0(i,j,1))
            tmp1 = squeeze(Qnet0(i,j,:))';
            % tmp1a = detrend(tmp1 - mean(tmp1));
            tmp1a = tmp1-mean(tmp1);
            
            %--------remove harmonic components ---------
            [p1,y1] = gft(t',tmp1a',kx,pvar,dvar,ifwd);
            tmp1a = tmp1a - y1';% remove harmonic signals
            Qnet_HighFreq(i,j,:)=tmp1a;
            Qnet_Harmonic(i,j,:)=y1;
            
            %-------Lowpass filtering the data---------------------
            L=length(tmp1a);
            tmp1a0(1:L)=flip(tmp1a);tmp1a0(L+1:2*L)=tmp1a;tmp1a0(2*L+1:3*L)=flip(tmp1a);
            f_width=19;
            tmp2a0=filter(hanning(f_width)/sum(hanning(f_width)),1,tmp1a0);
            tmp2b0=tmp2a0(L+(f_width-1)/2+1:2*L+(f_width-1)/2);
            Qnet_LowFreq(i,j,:)=tmp2b0;
        end
    end
end



save Qnet_LowFreq_Hanning  Qnet_LowFreq 
save Qnet_HighFreq_Hanning  Qnet_HighFreq 

%-----------Spatially smooth the data using box average--------------------
Lx=5;Ly=5; % need to adjust these scales for data have different spatial resolutions.  
for it=1:nt
   Qnet_LowFreq(:,:,it)=smooth2D_per(squeeze(Qnet_LowFreq(:,:,it)),Lx,Ly,lat);
end
save Qnet_LowFreq_Hanning_smooth  Qnet_LowFreq 

