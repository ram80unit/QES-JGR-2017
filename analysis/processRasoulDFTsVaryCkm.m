% read DFT results from Rasoul runs; compare ambient and perturbed.
% In this version, rather than comparing ambient vs. perturbed, we compare
% dipole (D) versus monopole (M) computed fields

clear all; close all;

topdir = '/Users/ram80unit/Google Drive/Projects/RasoulQES/510km_100Ckm/';
load([topdir 'alldfts.mat']);

h1 = figure(1);
ax(1) = subplot(411);
ax(2) = subplot(412);
ax(3) = subplot(413);
ax(4) = subplot(414);

set(0,'DefaultAxesFontSize',14);

nullwidth = 50; % grid cells to cut out near nulls
jumpthresh = 0.2;      % threshold of dphase/dx to consider a null
frind = 1;      % index of frequency to use; index 5 into Hp array is 24 kHz
mindist = 1300;     % distance range to analyze
maxdist = 2200;     % distance range to analyze

% data in alldfts.mat covers all of the runs shown as data_*.mat. Use the
% data file names as indices.

din = dir([topdir 'data_*.mat']);

return;

for m = 1:length(din),
    
    inputfile = din(m).name;
    
    i1 = find(dft(m).dist > mindist,1,'first');
    i2 = find(dft(m).dist > maxdist,1,'first');
    
    % ambient: these are the runs without any QES input.
    
    dist = dft(m).dist(1:i2);
    HpAmp0(:,m) = dft(m).Hp.amp(1:i2,frind);
    tempphase = dft(m).Hp.phase(1:i2,frind) * 180/pi;
    HpPhase0(:,m) = tempphase + (360/3e8 * dft(m).dftfreqs(frind)*dft(m).dist(1:i2)*1000);
    
    % perturbed
    
    HpAmp1(:,m) = dft(m+12).Hp.amp(1:i2,frind);
    tempphase = dft(m+12).Hp.phase(1:i2,frind) * 180/pi;
    HpPhase1(:,m) = tempphase + (360/3e8 * dft(m+12).dftfreqs(frind)*dft(m+12).dist(1:i2)*1000);
    
    % for the means, disregard 20 km regions past 1200 km where the diff of the
    % Ambient amplitude jumps by more than jumpthresh
    
   inds = find(abs(diff(20*log10(HpAmp1(:,m)))) > jumpthresh);
   for k = 1:length(inds),
       if inds(k) > nullwidth && inds(k) < (length(dist)-nullwidth),
           HpAmp1(inds(k)-nullwidth:inds(k)+nullwidth,m) = NaN;
           HpPhase1(inds(k)-nullwidth:inds(k)+nullwidth,m) = NaN;
       end
   end
    
end

% now, determine perturbations from HpAmp1 and HpPhase1 by comparing sets
% appropriately.

for m = 1:6,
    amppert(:,m) = 20*log10(HpAmp1(:,2*m-1)) - 20*log10(HpAmp1(:,2*m));
    phasepert(:,m) = HpPhase1(:,2*m-1) - HpPhase1(:,2*m);
end


inds = [1:6];
plot(ax(1),dist,20*log10(HpAmp1(:,inds)));

plot(ax(2),dist,amppert(:,inds));

plot(ax(3),dist,HpPhase1(:,inds));

plot(ax(4),dist,phasepert(:,inds));

xlabel(ax(4),'Distance (km)');
ylabel(ax(1),'Amplitude (dB)');
ylabel(ax(2),'Amplitude Perturbation (dB)');
ylabel(ax(3),'Phase (deg)');
ylabel(ax(4),'Phase perturbation (deg)');
title(ax(1),sprintf('Simulation %d of %d: %s',m,length(d)/2,inputfile));
hold(ax(1),'off');
hold(ax(3),'off');

set(ax(1),'xlim',[1200 max(dist)]);
set(ax(2),'xlim',[1200 max(dist)]);
set(ax(3),'xlim',[1200 max(dist)]);
set(ax(4),'xlim',[1200 max(dist)]);

% now, get mean amplitude and phase perturbation from chosen distance to end

for m = 1:6,
    
    inputfile = din(m).name;
    
    meanAmpPert(m) = mean(abs(amppert(i1:i2,m)),'omitnan');
    meanPhasePert(m) = mean(abs(phasepert(i1:i2,m)),'omitnan');
    ninetyAmp(m) = prctile(abs(amppert(i1:i2,m)),95);
    ninetyPhase(m) = prctile(abs(phasepert(i1:i2,m)),95);
    tenAmp(m) = prctile(abs(amppert(i1:i2,m)),5);
    tenPhase(m) = prctile(abs(phasepert(i1:i2,m)),5);
    
    fprintf('%s: DA = [%.3f %.3f %.3f] dB; Dphi = [%.3f %.3f %.3f] deg\n',inputfile,...
        tenAmp(m),meanAmpPert(m),ninetyAmp(m),tenPhase(m),meanPhasePert(m),ninetyPhase(m));
    
    
    pause(0.5);
    
end


%% now plot those results versus angle.

% there are five different ionospheres, and four different angles.

AmpPert.mean = meanAmpPert;
AmpPert.ninety = ninetyAmp;
AmpPert.ten = tenAmp;

PhasePert.mean = meanPhasePert;
PhasePert.ninety = ninetyPhase;
PhasePert.ten = tenPhase;


%% plot

h2 = figure(2);
set(h2,'position',[100 200 800 500]);
for m = 1:2,
    ax(m) = subplot(1,2,m);
end

ckm = [100 200 300 400 500 600];

err = [AmpPert.ninety fliplr(AmpPert.ten)];
errx = [ckm fliplr(ckm)];
patch(errx,err,'b','Parent',ax(1),'FaceAlpha',0.2,'EdgeAlpha',0);
hold(ax(1),'on');
plot(ax(1),ckm,AmpPert.mean,'b');
ylabel(ax(1),'Amp perturbation (dB)');

% plot phase perturbations too

err = [PhasePert.ninety fliplr(PhasePert.ten)];
errx = [ckm fliplr(ckm)];
patch(errx,err,'b','Parent',ax(2),'FaceAlpha',0.2,'EdgeAlpha',0);
hold(ax(2),'on');
plot(ax(2),ckm,PhasePert.mean,'b');
xlabel(ax(2),'Charge moment change (C-km)');
ylabel(ax(2),'Phase perturbation (deg)');
title(ax(2),'Ionosphere e1');
