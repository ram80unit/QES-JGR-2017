% plot Rasoul's perturbations

clear all; close all;

set(0,'DefaultAxesFontSize',16);

datadir = '/Users/Bob/Google Drive/Projects/RasoulQES/EarlyFastData1020/';

dD = dir([datadir '*_*D*.mat']);
dM = dir([datadir '*_*M*.mat']);

%inds = [11 12 13 15 3 8 13 23];
inds = [3 8 13 18 6 7 8 10];
%inds = [9:16];

h1 = figure(1);
set(h1,'position',[200 500 1500 700]);
for m = 1:8,
    ax(m) = subplot(2,4,m);
end

h2 = figure(2);
set(h2,'position',[300 400 1500 400]);
ax2(1) = subplot(1,3,1);
ax2(2) = subplot(1,3,2);
ax2(3) = subplot(1,3,3);

newx = [-150:1:150]*1e3;
newz = [0:1:120]*1e3;

for m = 1:length(inds),
    
    i = inds(m);
    
    A = load([datadir dD(i).name]);
    B = load([datadir dM(i).name]);
    
    [Nx,Nz] = meshgrid(A.x,A.z);
    nudipole = interp2(Nx,Nz,A.nu_eff',newx,newz');
    
    [Nx,Nz] = meshgrid(B.x,B.z);
    numonopole = interp2(Nx,Nz,B.nu_eff',newx,newz');
    
    bg = repmat(nudipole(:,1),1,size(nudipole,2));
    
    %ref = bg;
    ref = nudipole;
    
    diffnu = -(numonopole - ref)./ref * 100;
    
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
    caxis(ax(m),[-100 100]);
    set(ax(m),'xlim',[-120 120],'ylim',[0 100]);
    
    hold(ax(m),'on');
    plot(ax(m),[-120 120],[70 70],'k:');
    plot(ax(m),[-120 120],[75 75],'k:');
    plot(ax(m),[-120 120],[80 80],'k:');
    
    %pause;
    
end

set(ax2(1),'xlim',[-150 150],'ylim',[0 100]);
xlabel(ax2(1),'Distance (km)');
ylabel(ax2(1),'Change in nu (%)');
title(ax2(1),'Change in nu at 70 km');
legend(ax2(1),'Bdip = 10, e5','Bdip = 30, e5','Bdip = 60, e5','Bdip = 90, e5',...
    'Bdip = 30, e1','Bdip = 30, e3','Bdip = 30, e5','Bdip = 30, WS2');

set(ax2(2),'xlim',[-150 150],'ylim',[0 100]);
xlabel(ax2(2),'Distance (km)');
title(ax2(2),'Change in nu at 75 km');

set(ax2(3),'xlim',[-150 150],'ylim',[0 100]);
xlabel(ax2(3),'Distance (km)');
title(ax2(3),'Change in nu at 80 km');

colorbar('peer',ax(8),'East');

title(ax(1),'Bdip = 0, Ionosphere e5');
title(ax(2),'Bdip = 10, Ionosphere e5');
title(ax(3),'Bdip = 30, Ionosphere e5');
title(ax(4),'Bdip = 90, Ionosphere e5');

title(ax(5),'Bdip = 30, Ionosphere e1');
title(ax(6),'Bdip = 30, Ionosphere e3');
title(ax(7),'Bdip = 30, Ionosphere e5');
title(ax(8),'Bdip = 30, Ionosphere WS2');

xlabel(ax(5),'Distance (km)');
xlabel(ax(6),'Distance (km)');
xlabel(ax(7),'Distance (km)');
xlabel(ax(8),'Distance (km)');

ylabel(ax(1),'Altitude (km)');
ylabel(ax(5),'Altitude (km)');