% plot Rasoul's perturbations

clear all; close all;

set(0,'DefaultAxesFontSize',16);

%datadir = '/Users/Bob/Google Drive/Projects/RasoulQES/EarlyFastData510/';
datadir = '/Users/ram80unit/Google Drive/Projects/RasoulQES/510km_400Ckm/';

dM = dir([datadir '*_pCG.mat']);

inds = [1:8];

h1 = figure(1);
set(h1,'position',[200 500 1500 700]);
for m = 1:8,
    ax(m) = subplot(2,4,m);
end
c1 = colormap('jet');
set(0,'DefaultFigureColormap',c1);

h2 = figure(2);
set(h2,'position',[300 400 1500 400]);
ax2(1) = subplot(1,3,1);
ax2(2) = subplot(1,3,2);
ax2(3) = subplot(1,3,3);

newx = [-150:1:150]*1e3;
newz = [0:1:120]*1e3;

cmax = 30;  % maximum to show on colorscale

for m = 1:length(inds),
    
    i = inds(m);
    
    Mname = dM(i).name;
    Dname = [dM(i).name(1:end-8) '.mat'];
    Dname = strrep(Dname,'M','D');
    
    A = load([datadir Dname]);
    B = load([datadir Mname]);
    
    [Nx,Nz] = meshgrid(A.x,A.z);
    nudipole = interp2(Nx,Nz,A.nu_eff',newx,newz');
    
    [Nx,Nz] = meshgrid(B.x,B.z);
    numonopole = interp2(Nx,Nz,B.nu_eff',newx,newz');
    
    bg = repmat(nudipole(:,1),1,size(nudipole,2));
    
    %ref = bg;
    ref = nudipole;
    
    diffnu = (numonopole - ref)./ref * 100;
    
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
    
    hold(ax(m),'on');
    plot(ax(m),[-120 120],[70 70],'k:');
    plot(ax(m),[-120 120],[75 75],'k:');
    plot(ax(m),[-120 120],[80 80],'k:');
    
    legentry{m} = dM(i).name;
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

colorbar('peer',ax(8),'East');
colormap('jet');

for m = 1:8,
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
xlabel(ax(7),'Distance (km)');
xlabel(ax(8),'Distance (km)');

ylabel(ax(1),'Altitude (km)');
ylabel(ax(5),'Altitude (km)');