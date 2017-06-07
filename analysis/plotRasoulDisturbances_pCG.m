% plot Rasoul's perturbations

clear all; close all;

set(0,'DefaultAxesFontSize',16);

%datadir = '/Users/Bob/Google Drive/Projects/RasoulQES/EarlyFastData510/';
datadir = '/Users/Bob/Google Drive/Projects/RasoulQES/510km_400Ckm/';

dD = dir([datadir 'data_*D*.mat']);
dM = dir([datadir 'data*_pCG.mat']);

% need to carefully figure out indices to get pCG and backgrounds

inds = [3 8 13 18 23 11:15];

h1 = figure(1);
set(h1,'position',[200 500 1500 700]);
for m = 1:10,
    ax(m) = subplot(2,5,m);
end
c1 = colormap('jet');
set(0,'DefaultFigureColormap',c1);

h2 = figure(2);
set(h2,'position',[300 400 1500 700]);
for m = 1:10,
    ax2(m) = subplot(2,5,m);
end

newx = [-150:1:150]*1e3;
newz = [0:1:120]*1e3;

cmax = 80;  % maximum to show on colorscale

for m = 1:length(inds),
    
    i = inds(m);
    
    Mname = dM(i).name;
    Dname = dD(i).name;
    
    A = load([datadir Dname]);
    B = load([datadir Mname]);
    
    [Nx,Nz] = meshgrid(A.x,A.z);
    nudipole = interp2(Nx,Nz,A.nu_eff',newx,newz');
    
    [Nx,Nz] = meshgrid(B.x,B.z);
    numonopole = interp2(Nx,Nz,B.nu_eff',newx,newz');
    
    ref = nudipole;
    
    diffnu = (numonopole - ref)./ref * 100;
    
    altind = find(newz >= 70e3,1,'first');    
    diffnu1D70 = diffnu(altind,:);
    
    altind = find(newz >= 75e3,1,'first');    
    diffnu1D75 = diffnu(altind,:);
    
    altind = find(newz >= 80e3,1,'first');    
    diffnu1D80 = diffnu(altind,:);
    
    plot(ax2(m),newx/1e3,abs(diffnu1D70));
    hold(ax2(m),'on');
    
    plot(ax2(m),newx/1e3,abs(diffnu1D75));
    plot(ax2(m),newx/1e3,abs(diffnu1D80));
    
    set(ax2(m),'xlim',[-80 80],'ylim',[0 cmax]);
    title(ax2(m),'Change in nu at 70 km');
    legend(ax2(m),'70 km','75 km','80 km');
    
    imagesc(newx/1e3,newz/1e3,diffnu,'parent',ax(m)); axis(ax(m),'xy');
    caxis(ax(m),[-cmax cmax]);
    set(ax(m),'xlim',[-120 120],'ylim',[0 100]);
    
    hold(ax(m),'on');
    plot(ax(m),[-120 120],[70 70],'k:');
    plot(ax(m),[-120 120],[75 75],'k:');
    plot(ax(m),[-120 120],[80 80],'k:');
    
    legentry{m} = dM(i).name;
    title(ax(m),legentry{m});
    title(ax2(m),legentry{m});
    %pause;
    
end

colorbar('peer',ax(10),'East');
colormap('jet');

%title(ax(1),'Bdip = 10, Ionosphere e5');
%title(ax(2),'Bdip = 30, Ionosphere e5');
%title(ax(3),'Bdip = 60, Ionosphere e5');
%title(ax(4),'Bdip = 90, Ionosphere e5');

%title(ax(5),'Bdip = 30, Ionosphere e1');
%title(ax(6),'Bdip = 30, Ionosphere e3');
%title(ax(7),'Bdip = 30, Ionosphere e5');
%title(ax(8),'Bdip = 30, Ionosphere WS2');

xlabel(ax(6),'Distance (km)');
xlabel(ax(7),'Distance (km)');
xlabel(ax(8),'Distance (km)');
xlabel(ax(9),'Distance (km)');
xlabel(ax(10),'Distance (km)');

ylabel(ax(1),'Altitude (km)');
ylabel(ax(6),'Altitude (km)');

xlabel(ax2(6),'Distance (km)');
xlabel(ax2(7),'Distance (km)');
xlabel(ax2(8),'Distance (km)');
xlabel(ax2(9),'Distance (km)');
xlabel(ax2(10),'Distance (km)');

ylabel(ax2(1),'Change in nu (%)');
ylabel(ax2(6),'Change in nu (%)');