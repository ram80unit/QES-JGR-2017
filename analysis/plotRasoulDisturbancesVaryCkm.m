% plot Rasoul's perturbations

clear all; close all;

set(0,'DefaultAxesFontSize',16);

datadir = '/Users/ram80unit/Google Drive/Projects/RasoulQES/alt_var_const_charge_moment/';

dD = dir([datadir 'data_30D*.mat']);
dM = dir([datadir 'data_30M*.mat']);


h1 = figure(1);
set(h1,'position',[200 500 1500 700]);
for m = 1:6,
    ax(m) = subplot(2,3,m);
end
c1 = colormap('jet');
set(0,'DefaultFigureColormap',c1);

h2 = figure(2);
set(h2,'position',[300 400 1500 400]);
ax2(1) = subplot(1,3,1);
ax2(2) = subplot(1,3,2);
ax2(3) = subplot(1,3,3);

h3 = figure(3);
set(h3,'position',[200 500 1500 700]);
for m = 1:6,
    ax3(m) = subplot(2,3,m);
end


newx = [-150:1:150]*1e3;
newz = [0:1:120]*1e3;

cmax = 80;  % maximum to show on colorscale

for m = 1:6,
    
    i = m;
    
    A = load([datadir dD(i).name]);
    B = load([datadir dM(i).name]);
    
    [Nx,Nz] = meshgrid(A.x,A.z);
    nudipole = interp2(Nx,Nz,A.nu_eff',newx,newz');
    nedipole = interp2(Nx,Nz,A.ne',newx,newz');
    
    [Nx,Nz] = meshgrid(B.x,B.z);
    numonopole = interp2(Nx,Nz,B.nu_eff',newx,newz');
    nemonopole = interp2(Nx,Nz,B.ne',newx,newz');

    bg = repmat(nudipole(:,1),1,size(nudipole,2));
    bgne = repmat(nedipole(:,1),1,size(nedipole,2));
    
    %ref = bg;
    ref = nudipole;
    refne = nedipole;
    
    diffnu = (numonopole - ref)./ref * 100;
    diffne = (nemonopole - refne)./refne * 100;
    
    altind = find(newz >= 70e3,1,'first');    
    diffnu1D70 = diffnu(altind,:);
    
    altind = find(newz >= 75e3,1,'first');    
    diffnu1D75 = diffnu(altind,:);
    
    altind = find(newz >= 80e3,1,'first');    
    diffnu1D80 = diffnu(altind,:);
    
    plot(ax2(1),newx/1e3,abs(diffnu1D70));
    hold(ax2(1),'on');
    
    plot(ax2(2),newx/1e3,abs(diffnu1D75));
    hold(ax2(2),'on');
    
    plot(ax2(3),newx/1e3,abs(diffnu1D80));
    hold(ax2(3),'on');
    
    imagesc(newx/1e3,newz/1e3,diffnu,'parent',ax(m)); axis(ax(m),'xy');
    caxis(ax(m),[-cmax cmax]);
    set(ax(m),'xlim',[-120 120],'ylim',[0 100]);
    
    imagesc(newx/1e3,newz/1e3,diffne,'parent',ax3(m)); axis(ax3(m),'xy');
    caxis(ax3(m),[-cmax cmax]);
    set(ax3(m),'xlim',[-120 120],'ylim',[0 100]);
    
    hold(ax(m),'on');
    plot(ax(m),[-120 120],[70 70],'k:');
    plot(ax(m),[-120 120],[75 75],'k:');
    plot(ax(m),[-120 120],[80 80],'k:');
    
    legentry{m} = dD(i).name;
    %pause;
    
end

set(ax2(1),'xlim',[-80 80],'ylim',[0 cmax]);
xlabel(ax2(1),'Distance (km)');
ylabel(ax2(1),'Change in nu (%)');
title(ax2(1),'Change in nu at 70 km');
legend(ax2(1),legentry);

set(ax2(2),'xlim',[-80 80],'ylim',[0 cmax]);
xlabel(ax2(2),'Distance (km)');
title(ax2(2),'Change in nu at 75 km');

set(ax2(3),'xlim',[-80 80],'ylim',[0 cmax]);
xlabel(ax2(3),'Distance (km)');
title(ax2(3),'Change in nu at 80 km');

colorbar('peer',ax(6),'East');
colormap('jet');

for m = 1:6,
    title(ax(m),legentry{m});
end

%title(ax(1),'Bdip = 10, Ionosphere e5');
%title(ax(2),'Bdip = 30, Ionosphere e5');
%title(ax(3),'Bdip = 60, Ionosphere e5');
%title(ax(4),'Bdip = 90, Ionosphere e5');

%title(ax(5),'Bdip = 30, Ionosphere e1');
%title(ax(6),'Bdip = 30, Ionosphere e3');
%title(ax(7),'Bdip = 30, Ionosphere e5');
%title(ax(8),'Bdip = 30, Ionosphere WS2');

xlabel(ax(5),'Distance (km)');
xlabel(ax(6),'Distance (km)');

ylabel(ax(1),'Altitude (km)');
ylabel(ax(5),'Altitude (km)');