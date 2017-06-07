% read DFT results from Rasoul runs; compare ambient and perturbed.
% In this version, rather than comparing ambient vs. perturbed, we compare
% dipole (D) versus monopole (M) computed fields

clear all; close all;

topdir = '/Users/Bob/Google Drive/Projects/RasoulQES/NewAltVariation/';
load([topdir 'alldfts.mat']);

nullwidth = 50; % grid cells to cut out near nulls
jumpthresh = 0.2;      % threshold of dphase/dx to consider a null
frind = 3;      % index of frequency to use; index 5 into Hp array is 24 kHz
mindist = 1300;     % distance range to analyze
maxdist = 2200;     % distance range to analyze

% indices of distance of interest    
i1 = find(dft(1).dist > mindist,1,'first');
i2 = find(dft(1).dist > maxdist,1,'first');
    
% is these runs, indices 1-5 are for a lower charge at 10 km, upper at
% 12:2:20.
% indices 6-8 are for lower at 15 km, upper at 16:2:20.
% indices 9:13 are for lower at 5 km, upper at 12:2:20.
% indices 14:16 are the monopoles for those sets.

% let's pull them out one at a time.


for m = 1:7,
    
    % ambient
    
    dist = dft(m).dist(1:i2);
    HpAmp0(:,m) = dft(m).Hp.amp(1:i2,frind);
    tempphase = dft(m).Hp.phase(1:i2,frind) * 180/pi;
    HpPhase0(:,m) = tempphase + (360/3e8 * dft(m).dftfreqs(frind)*dft(m).dist(1:i2)*1000);
    
    % perturbed
    
    dist = dft(m+7).dist(1:i2);
    HpAmp1(:,m) = dft(m+7).Hp.amp(1:i2,frind);
    tempphase = dft(m+7).Hp.phase(1:i2,frind) * 180/pi;
    HpPhase1(:,m) = tempphase + (360/3e8 * dft(m+7).dftfreqs(frind)*dft(m+7).dist(1:i2)*1000);
    
    % for the means, disregard 20 km regions past 1200 km where the diff of the
    % Ambient amplitude jumps by more than jumpthresh
    
    inds = find(abs(diff(20*log10(HpAmp1(:,m)))) > jumpthresh);
    for k = 1:length(inds),
        if inds(k) > nullwidth && inds(k) < (length(dist)-nullwidth),
            HpAmp1(inds(k)-nullwidth:inds(k)+nullwidth,m) = NaN;
            HpPhase1(inds(k)-nullwidth:inds(k)+nullwidth,m) = NaN;
        end
    end
    
    % perturbations
    amppert(:,m) = 20*log10(HpAmp1(:,m)) - 20*log10(HpAmp0(:,m));
    phasepert(:,m) = HpPhase1(:,m) - HpPhase0(:,m);
    
end

%% now, get mean amplitude and phase perturbation from chosen distance to end

for m = 1:7,
    
    
    AmpPert50(m) = prctile(abs(amppert(i1:i2,m)),50);
    PhasePert50(m) = prctile(abs(phasepert(i1:i2,m)),50);
    AmpPert90(m) = prctile(abs(amppert(i1:i2,m)),90);
    PhasePert90(m) = prctile(abs(phasepert(i1:i2,m)),90);
    AmpPert10(m) = prctile(abs(amppert(i1:i2,m)),10);
    PhasePert10(m) = prctile(abs(phasepert(i1:i2,m)),10);
    
    fprintf('%d: DA = [%.3f %.3f %.3f] dB; Dphi = [%.3f %.3f %.3f] deg\n',m,...
        AmpPert10(m),AmpPert50(m),AmpPert90(m),PhasePert10(m),PhasePert50(m),PhasePert90(m));
    
    
    pause(0.5);
    
end

%% 

h1 = figure(1);
ax(1) = subplot(411);
ax(2) = subplot(412);
ax(3) = subplot(413);
ax(4) = subplot(414);

set(0,'DefaultAxesFontSize',14);

inds = [1:7];
plot(ax(1),dist,20*log10(HpAmp1(:,inds)));

plot(ax(2),dist,amppert(:,inds));

plot(ax(3),dist,HpPhase1(:,inds));

plot(ax(4),dist,phasepert(:,inds));

xlabel(ax(4),'Distance (km)');
ylabel(ax(1),'Amplitude (dB)');
ylabel(ax(2),'Amplitude Perturbation (dB)');
ylabel(ax(3),'Phase (deg)');
ylabel(ax(4),'Phase perturbation (deg)');
%title(ax(1),sprintf('Simulation %d of %d: %s',m,length(d)/2,inputfile));
hold(ax(1),'off');
hold(ax(3),'off');

set(ax(1),'xlim',[1200 max(dist)]);
set(ax(2),'xlim',[1200 max(dist)]);
set(ax(3),'xlim',[1200 max(dist)]);
set(ax(4),'xlim',[1200 max(dist)]);

%% now plot those results versus upper charge layer height.

% there are five different ionospheres, and four different angles.

AmpPert.mean = AmpPert50([7 1:6]);
AmpPert.ninety = AmpPert90([7 1:6]);
AmpPert.ten = AmpPert10([7 1:6]);

% phase

PhasePert.mean = PhasePert50([7 1:6]);
PhasePert.ninety = PhasePert90([7 1:6]);
PhasePert.ten = PhasePert10([7 1:6]);



%% plot

h2 = figure(2);
set(h2,'position',[100 495 1600 400]);
ax(1) = subplot(1,2,1);
ax(2) = subplot(1,2,2);

angle = [8:2:20];

err = [AmpPert.ninety fliplr(AmpPert.ten)];
errx = [angle fliplr(angle)];
patch(errx,err,'b','Parent',ax(1),'FaceAlpha',0.2,'EdgeAlpha',0);
hold(ax(1),'on');
plot(ax(1),angle,AmpPert.mean,'b');
set(ax(1),'xlim',[12 20]);
xlabel(ax(1),'Upper layer altitude (km)');
ylabel(ax(1),'Amp perturbation (dB)');
title(ax(1),'Lower layer at 5 km; constant Ckm');

err = [PhasePert.ninety fliplr(PhasePert.ten)];
errx = [angle fliplr(angle)];
patch(errx,err,'b','Parent',ax(2),'FaceAlpha',0.2,'EdgeAlpha',0);
hold(ax(2),'on');
plot(ax(2),angle,PhasePert.mean,'b');
set(ax(2),'xlim',[12 20]);
xlabel(ax(2),'Upper layer altitude (km)');
ylabel(ax(2),'Phase perturbation (deg)');
title(ax(2),'Lower layer at 5 km; constant Ckm');

set(ax(1),'ylim',[0 0.5]);
set(ax(2),'ylim',[0 3]);
