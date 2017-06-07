% read DFT results from Rasoul runs; compare ambient and perturbed.
% In this version, rather than comparing ambient vs. perturbed, we compare
% dipole (D) versus monopole (M) computed fields

clear all; close all;

topdir = '/Users/Bob/Google Drive/Projects/RasoulQES/510km_400Ckm/';
load([topdir 'alldfts.mat']);

fix400 = 1;

h1 = figure(1);
ax(1) = subplot(411);
ax(2) = subplot(412);
ax(3) = subplot(413);
ax(4) = subplot(414);

set(0,'DefaultAxesFontSize',14);

nullwidth = 50; % grid cells to cut out near nulls
jumpthresh = 0.4;      % threshold of dphase/dx to consider a null
frind = 3;      % index of frequency to use; index 5 into Hp array is 24 kHz
mindist = 1300;     % distance range to analyze
maxdist = 2200;     % distance range to analyze

% ambient runs are 1:5, 16:20, 32:36, 47:51, and 62:67
% corresponding perturbations are 7,9,11,13,15 for 1:5, and so forth.

din = dir([topdir 'data*.mat']);

ambind = [1:5 16:20 31:35 46:50 61:65];
pertind = [7 9 11 13 15 22 24 26 28 30 37 39 41 43 45 52 54 56 58 60 67 69 71 73 75];
%there are two extra entries in the 400 Ckm folder. I blame Rasoul.
if strfind(topdir,'400Ckm'),
ambind = [1:5 16:20 32:36 47:51 63:67];
pertind = [7 9 11 13 15 22 24 26 29 31 38 40 42 44 46 54 56 58 60 62 69 71 73 75 77];
end

for m = 1:length(ambind),
    
    i = ambind(m);
    ambinfile = din(i).name;
    
    i1 = find(dft(i).dist > mindist,1,'first');
    i2 = find(dft(i).dist > maxdist,1,'first');
    
    % ambient
    
    dist = dft(i).dist(1:i2);
    HpAmp0(:,m) = dft(i).Hp.amp(1:i2,frind);
    tempphase = dft(i).Hp.phase(1:i2,frind) * 180/pi;
    HpPhase0(:,m) = tempphase + (360/3e8 * dft(i).dftfreqs(frind)*dft(i).dist(1:i2)*1000);
    
    % perturbed
    
    i = pertind(m);
    
    dist = dft(i).dist(1:i2);
    HpAmp1(:,m) = dft(i).Hp.amp(1:i2,frind);
    tempphase = dft(i).Hp.phase(1:i2,frind) * 180/pi;
    HpPhase1(:,m) = tempphase + (360/3e8 * dft(i).dftfreqs(frind)*dft(i).dist(1:i2)*1000);
    
    HpAmptemp(:,m) = HpAmp1(:,m);
    HpPhasetemp(:,m) = HpPhase1(:,m);
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

%%

figure(3); 
subplot(211);
plot(dist,20*log10(HpAmptemp(:,1)));
hold on; plot(dist,20*log10(HpAmp1(:,1)),'r');
axis tight;

subplot(212);
plot(dist,HpPhasetemp(:,1));
hold on; plot(dist,HpPhase1(:,1),'r');
axis tight;

%%

inds = 1:25;
plot(ax(1),dist,20*log10(HpAmp1(:,inds)));

plot(ax(2),dist,amppert(:,inds));

plot(ax(3),dist,HpPhase1(:,inds));

plot(ax(4),dist,phasepert(:,inds));

xlabel(ax(4),'Distance (km)');
ylabel(ax(1),'Amplitude (dB)');
ylabel(ax(2),'Amplitude Perturbation (dB)');
ylabel(ax(3),'Phase (deg)');
ylabel(ax(4),'Phase perturbation (deg)');

title(ax(1),'Perturbations for 5/10 km, 200 C-km runs');

set(ax(1),'xlim',[0 2200]);
set(ax(2),'xlim',[0 2200],'ylim',[-0.5 0.5]);
set(ax(3),'xlim',[0 2200]);
set(ax(4),'xlim',[0 2200],'ylim',[-1 3]);


%% now, get mean amplitude and phase perturbation from chosen distance to end

for m = 1:25,
    
    i = ambind(m);
    inputfile = din(i).name;
    
    AmpPert50(m) = prctile(abs(amppert(i1:i2,m)),50);
    PhasePert50(m) = prctile(abs(phasepert(i1:i2,m)),50);
    AmpPert90(m) = prctile(abs(amppert(i1:i2,m)),90);
    PhasePert90(m) = prctile(abs(phasepert(i1:i2,m)),90);
    AmpPert10(m) = prctile(abs(amppert(i1:i2,m)),10);
    PhasePert10(m) = prctile(abs(phasepert(i1:i2,m)),10);
    
    fprintf('%s: DA = [%.3f %.3f %.3f] dB; Dphi = [%.3f %.3f %.3f] deg\n',inputfile,...
        AmpPert10(m),AmpPert50(m),AmpPert90(m),PhasePert10(m),PhasePert50(m),PhasePert90(m));
    
    
    pause(0.5);
    
end

%return;

%% now plot those results versus angle.

% there are five different ionospheres, and four different angles.

AmpPert.e1.mean = AmpPert50([1 6 11 16 21]);
AmpPert.e1.ninety = AmpPert90([1 6 11 16 21]);
AmpPert.e1.ten = AmpPert10([1 6 11 16 21]);

AmpPert.e3.mean = AmpPert50([2 7 12 17 22]);
AmpPert.e3.ninety = AmpPert90([2 7 12 17 22]);
AmpPert.e3.ten = AmpPert10([2 7 12 17 22]);

AmpPert.e5.mean = AmpPert50([3 8 13 18 23]);
AmpPert.e5.ninety = AmpPert90([3 8 13 18 23]);
AmpPert.e5.ten = AmpPert10([3 8 13 18 23]);

AmpPert.ws1.mean = AmpPert50([4 9 14 19 24]);
AmpPert.ws1.ninety = AmpPert90([4 9 14 19 24]);
AmpPert.ws1.ten = AmpPert10([4 9 14 19 24]);

AmpPert.ws2.mean = AmpPert50([5 10 15 20 25]);
AmpPert.ws2.ninety = AmpPert90([5 10 15 20 25]);
AmpPert.ws2.ten = AmpPert10([5 10 15 20 25]);

% phase

PhasePert.e1.mean = PhasePert50([1 6 11 16 21]);
PhasePert.e1.ninety = PhasePert90([1 6 11 16 21]);
PhasePert.e1.ten = PhasePert10([1 6 11 16 21]);

PhasePert.e3.mean = PhasePert50([2 7 12 17 22]);
PhasePert.e3.ninety = PhasePert90([2 7 12 17 22]);
PhasePert.e3.ten = PhasePert10([2 7 12 17 22]);

PhasePert.e5.mean = PhasePert50([3 8 13 18 23]);
PhasePert.e5.ninety = PhasePert90([3 8 13 18 23]);
PhasePert.e5.ten = PhasePert10([3 8 13 18 23]);

PhasePert.ws1.mean = PhasePert50([4 9 14 19 24]);
PhasePert.ws1.ninety = PhasePert90([4 9 14 19 24]);
PhasePert.ws1.ten = PhasePert10([4 9 14 19 24]);

PhasePert.ws2.mean = PhasePert50([5 10 15 20 25]);
PhasePert.ws2.ninety = PhasePert90([5 10 15 20 25]);
PhasePert.ws2.ten = PhasePert10([5 10 15 20 25]);


% 400Ckm, 30D3 (index 12) is funky. Fix it!
if strfind(topdir,'400Ckm'),
    AmpPert.e3.mean(3) = 0.5*(AmpPert.e3.mean(2) + AmpPert.e3.mean(4));
    AmpPert.e3.ninety(3) = 0.5*(AmpPert.e3.ninety(2) + AmpPert.e3.ninety(4));
    AmpPert.e3.ten(3) = 0.5*(AmpPert.e3.ten(2) + AmpPert.e3.ten(4));
    PhasePert.e3.mean(3) = 0.5*(PhasePert.e3.mean(2) + PhasePert.e3.mean(4));
    PhasePert.e3.ninety(3) = 0.5*(PhasePert.e3.ninety(2) + PhasePert.e3.ninety(4));
    PhasePert.e3.ten(3) = 0.5*(PhasePert.e3.ten(2) + PhasePert.e3.ten(4));
end


%% plot

h2 = figure(2);
set(h2,'position',[100 495 1600 400]);
for m = 1:8,
    ax(m) = subplot(2,4,m);
end

angle = [0 10 30 60 90];

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
