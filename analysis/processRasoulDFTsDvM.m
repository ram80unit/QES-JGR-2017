% read DFT results from Rasoul runs; compare ambient and perturbed.
% In this version, rather than comparing ambient vs. perturbed, we compare
% dipole (D) versus monopole (M) computed fields

clear all; close all;

topdir = '/Users/Bob/Google Drive/Projects/RasoulQES/510km_100Ckm/';
load([topdir 'alldfts.mat']);

h1 = figure(1);
ax(1) = subplot(411);
ax(2) = subplot(412);
ax(3) = subplot(413);
ax(4) = subplot(414);

set(0,'DefaultAxesFontSize',14);

nullwidth = 50; % grid cells to cut out near nulls
jumpthresh = 0.2;      % threshold of dphase/dx to consider a null
frind = 5;      % index of frequency to use; index 5 into Hp array is 24 kHz
mindist = 1300;     % distance range to analyze
maxdist = 2200;     % distance range to analyze

% first forty runs are ambient; next forty are perturbed. So, need to read
% into dft array starting at 40, and then offsetting appropriately.

din = dir([topdir 'data*.mat']);

for m = 1:length(din),
    
    inputfile = din(m).name;
    
    i1 = find(dft(m).dist > mindist,1,'first');
    i2 = find(dft(m).dist > maxdist,1,'first');
    
    % ambient
    
    dist = dft(m).dist(1:i2);
    HpAmp0(:,m) = dft(m).Hp.amp(1:i2,frind);
    tempphase = dft(m).Hp.phase(1:i2,frind) * 180/pi;
    HpPhase0(:,m) = tempphase + (360/3e8 * dft(m).dftfreqs(frind)*dft(m).dist(1:i2)*1000);
    
    % perturbed
    
    HpAmp1(:,m) = dft(m+40).Hp.amp(1:i2,frind);
    tempphase = dft(m+40).Hp.phase(1:i2,frind) * 180/pi;
    HpPhase1(:,m) = tempphase + (360/3e8 * dft(m+40).dftfreqs(frind)*dft(m+40).dist(1:i2)*1000);
    
    % for the means, disregard 20 km regions past 1200 km where the diff of the
    % Ambient amplitude jumps by more than jumpthresh
    
%     inds = find(abs(diff(20*log10(HpAmp1(:,m)))) > jumpthresh);
%     for k = 1:length(inds),
%         if inds(k) > nullwidth && inds(k) < (length(dist)-nullwidth),
%             HpAmp1(inds(k)-nullwidth:inds(k)+nullwidth,m) = NaN;
%             HpPhase1(inds(k)-nullwidth:inds(k)+nullwidth,m) = NaN;
%         end
%     end
    
end

% now, determine perturbations from HpAmp1 and HpPhase1 by comparing sets
% appropriately.

for m = 1:5,
    amppert(:,m) = 20*log10(HpAmp1(:,m)) - 20*log10(HpAmp1(:,m+5));
    phasepert(:,m) = HpPhase1(:,m) - HpPhase1(:,m+5);
end
for m = 6:10,
    amppert(:,m) = 20*log10(HpAmp1(:,m+5)) - 20*log10(HpAmp1(:,m+10));
    phasepert(:,m) = HpPhase1(:,m+5) - HpPhase1(:,m+10);
end
for m = 11:15,
    amppert(:,m) = 20*log10(HpAmp1(:,m+10)) - 20*log10(HpAmp1(:,m+15));
    phasepert(:,m) = HpPhase1(:,m+10) - HpPhase1(:,m+15);
end
for m = 16:20,
    amppert(:,m) = 20*log10(HpAmp1(:,m+15)) - 20*log10(HpAmp1(:,m+20));
    phasepert(:,m) = HpPhase1(:,m+15) - HpPhase1(:,m+20);
end


plot(ax(1),dist,20*log10(HpAmp1));

plot(ax(2),dist,amppert);

plot(ax(3),dist,HpPhase1);

plot(ax(4),dist,phasepert);

xlabel(ax(4),'Distance (km)');
ylabel(ax(1),'Amplitude (dB)');
ylabel(ax(2),'Amplitude Perturbation (dB)');
ylabel(ax(3),'Phase (deg)');
ylabel(ax(4),'Phase perturbation (deg)');
title(ax(1),sprintf('Simulation %d of %d: %s',m,length(d)/2,inputfile));
hold(ax(1),'off');
hold(ax(3),'off');

set(ax(1),'xlim',[0 max(dist)]);
set(ax(2),'xlim',[0 max(dist)]);
set(ax(3),'xlim',[0 max(dist)]);
set(ax(4),'xlim',[0 max(dist)]);

% now, get mean amplitude and phase perturbation from chosen distance to end

for m = 1:20,
    
    inputfile = din(m).name;
    
    meanAmpPert(m) = mean(abs(amppert(i1:i2,m)),'omitnan');
    meanPhasePert(m) = mean(abs(phasepert(i1:i2,m)),'omitnan');
    ninetyAmp(m) = prctile(abs(amppert(i1:i2,m)),90);
    ninetyPhase(m) = prctile(abs(phasepert(i1:i2,m)),90);
    tenAmp(m) = prctile(abs(amppert(i1:i2,m)),10);
    tenPhase(m) = prctile(abs(phasepert(i1:i2,m)),10);
    
    fprintf('%s: DA = [%.3f %.3f %.3f] dB; Dphi = [%.3f %.3f %.3f] deg\n',inputfile,...
        tenAmp(m),meanAmpPert(m),ninetyAmp(m),tenPhase(m),meanPhasePert(m),ninetyPhase(m));
    
    
    pause(0.5);
    
end


%% now plot those results versus angle.

% there are five different ionospheres, and four different angles.

AmpPert.e1.mean = meanAmpPert([1 6 11 16]);
AmpPert.e1.ninety = ninetyAmp([1 6 11 16]);
AmpPert.e1.ten = tenAmp([1 6 11 16]);

AmpPert.e3.mean = meanAmpPert([2 7 12 17]);
AmpPert.e3.ninety = ninetyAmp([2 7 12 17]);
AmpPert.e3.ten = tenAmp([2 7 12 17]);

AmpPert.e5.mean = meanAmpPert([3 8 13 18]);
AmpPert.e5.ninety = ninetyAmp([3 8 13 18]);
AmpPert.e5.ten = tenAmp([3 8 13 18]);

AmpPert.ws1.mean = meanAmpPert([4 9 14 19]);
AmpPert.ws1.ninety = ninetyAmp([4 9 14 19]);
AmpPert.ws1.ten = tenAmp([4 9 14 19]);

AmpPert.ws2.mean = meanAmpPert([5 10 15 20]);
AmpPert.ws2.ninety = ninetyAmp([5 10 15 20]);
AmpPert.ws2.ten = tenAmp([5 10 15 20]);

% phase

PhasePert.e1.mean = meanPhasePert([1 6 11 16]);
PhasePert.e1.ninety = ninetyPhase([1 6 11 16]);
PhasePert.e1.ten = tenPhase([1 6 11 16]);

PhasePert.e3.mean = meanPhasePert([2 7 12 17]);
PhasePert.e3.ninety = ninetyPhase([2 7 12 17]);
PhasePert.e3.ten = tenPhase([2 7 12 17]);

PhasePert.e5.mean = meanPhasePert([3 8 13 18]);
PhasePert.e5.ninety = ninetyPhase([3 8 13 18]);
PhasePert.e5.ten = tenPhase([3 8 13 18]);

PhasePert.ws1.mean = meanPhasePert([4 9 14 19]);
PhasePert.ws1.ninety = ninetyPhase([4 9 14 19]);
PhasePert.ws1.ten = tenPhase([4 9 14 19]);

PhasePert.ws2.mean = meanPhasePert([5 10 15 20]);
PhasePert.ws2.ninety = ninetyPhase([5 10 15 20]);
PhasePert.ws2.ten = tenPhase([5 10 15 20]);

%% plot

h2 = figure(2);
set(h2,'position',[100 495 1600 400]);
for m = 1:8,
    ax(m) = subplot(2,4,m);
end

angle = [10 30 60 90];

err = [AmpPert.e1.ninety fliplr(AmpPert.e1.ten)];
errx = [angle fliplr(angle)];
patch(errx,err,'b','Parent',ax(1),'FaceAlpha',0.2,'EdgeAlpha',0);
hold(ax(1),'on');
plot(ax(1),angle,AmpPert.e1.mean,'b');
set(ax(1),'xlim',[0 90]);
ylabel(ax(1),'Amp perturbation (dB)');
title(ax(1),'Ionosphere e1');

err = [AmpPert.e3.ninety fliplr(AmpPert.e3.ten)];
errx = [angle fliplr(angle)];
patch(errx,err,'b','Parent',ax(2),'FaceAlpha',0.2,'EdgeAlpha',0);
hold(ax(2),'on');
plot(ax(2),angle,AmpPert.e3.mean,'b');
set(ax(2),'xlim',[0 90]);
title(ax(2),'Ionosphere e3');

err = [AmpPert.e5.ninety fliplr(AmpPert.e5.ten)];
errx = [angle fliplr(angle)];
patch(errx,err,'b','Parent',ax(3),'FaceAlpha',0.2,'EdgeAlpha',0);
hold(ax(3),'on');
plot(ax(3),angle,AmpPert.e5.mean,'b');
set(ax(3),'xlim',[0 90]);
title(ax(3),'Ionosphere e5');

% err = [AmpPert.ws1.ninety fliplr(AmpPert.ws1.ten)];
% errx = [angle fliplr(angle)];
% patch(errx,err,'b','Parent',ax(4),'FaceAlpha',0.2,'EdgeAlpha',0);
% hold(ax(4),'on');
% plot(ax(4),angle,AmpPert.ws1.mean,'b');
% set(ax(4),'xlim',[0 90]);
% title(ax(4),'Ionosphere WS1');

err = [AmpPert.ws2.ninety fliplr(AmpPert.ws2.ten)];
errx = [angle fliplr(angle)];
patch(errx,err,'b','Parent',ax(4),'FaceAlpha',0.2,'EdgeAlpha',0);
hold(ax(4),'on');
plot(ax(4),angle,AmpPert.ws2.mean,'b');
set(ax(4),'xlim',[0 90]);
title(ax(4),'Ionosphere WS2');

% plot phase perturbations too

err = [PhasePert.e1.ninety fliplr(PhasePert.e1.ten)];
errx = [angle fliplr(angle)];
patch(errx,err,'b','Parent',ax(5),'FaceAlpha',0.2,'EdgeAlpha',0);
hold(ax(5),'on');
plot(ax(5),angle,PhasePert.e1.mean,'b');
set(ax(5),'xlim',[0 90]);
xlabel(ax(5),'B dip angle (deg)');
ylabel(ax(5),'Phase perturbation (deg)');
title(ax(5),'Ionosphere e1');

err = [PhasePert.e3.ninety fliplr(PhasePert.e3.ten)];
errx = [angle fliplr(angle)];
patch(errx,err,'b','Parent',ax(6),'FaceAlpha',0.2,'EdgeAlpha',0);
hold(ax(6),'on');
plot(ax(6),angle,PhasePert.e3.mean,'b');
set(ax(6),'xlim',[0 90]);
xlabel(ax(6),'B dip angle (deg)');
title(ax(6),'Ionosphere e3');

err = [PhasePert.e5.ninety fliplr(PhasePert.e5.ten)];
errx = [angle fliplr(angle)];
patch(errx,err,'b','Parent',ax(7),'FaceAlpha',0.2,'EdgeAlpha',0);
hold(ax(7),'on');
plot(ax(7),angle,PhasePert.e5.mean,'b');
set(ax(7),'xlim',[0 90]);
xlabel(ax(7),'B dip angle (deg)');
title(ax(7),'Ionosphere e5');

% err = [PhasePert.ws1.ninety fliplr(PhasePert.ws1.ten)];
% errx = [angle fliplr(angle)];
% patch(errx,err,'b','Parent',ax(9),'FaceAlpha',0.2,'EdgeAlpha',0);
% hold(ax(9),'on');
% plot(ax(9),angle,PhasePert.ws1.mean,'b');
% set(ax(9),'xlim',[0 90]);
% xlabel(ax(9),'B dip angle (deg)');
% title(ax(9),'Ionosphere WS1');

err = [PhasePert.ws2.ninety fliplr(PhasePert.ws2.ten)];
errx = [angle fliplr(angle)];
patch(errx,err,'b','Parent',ax(8),'FaceAlpha',0.2,'EdgeAlpha',0);
hold(ax(8),'on');
plot(ax(8),angle,PhasePert.ws2.mean,'b');
set(ax(8),'xlim',[0 90]);
xlabel(ax(8),'B dip angle (deg)');
title(ax(8),'Ionosphere WS2');