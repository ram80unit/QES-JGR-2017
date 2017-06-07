% plot DFT outputs from EMP code

function out = plotDFTemp(rundir,doplot)

   %doplot = 1;

%rundir = '/Users/ram80unit/Google Drive/STOIC/runs/emp/I0_5e04/';

s = get2drunparams(rundir,'double');

fid = fopen([rundir 'dft.dat'],'r');
ndft = fread(fid,1,'int');
dftfreqs = fread(fid,ndft,'double');
dft.Er = fread(fid,[2*ndft s.hh],'double');
dft.Et = fread(fid,[2*ndft s.hh],'double');
dft.Ep = fread(fid,[2*ndft s.hh],'double');
dft.Hr = fread(fid,[2*ndft s.hh],'double');
dft.Ht = fread(fid,[2*ndft s.hh],'double');
dft.Hp = fread(fid,[2*ndft s.hh],'double');
fclose(fid);

out.Er.amp = zeros(s.hh,ndft);
out.Er.phase = zeros(s.hh,ndft);
%out.Et.amp = zeros(s.hh,ndft);
%out.Et.phase = zeros(s.hh,ndft);
%out.Ep.amp = zeros(s.hh,ndft);
%out.Ep.phase = zeros(s.hh,ndft);
%out.Hr.amp = zeros(s.hh,ndft);
%out.Hr.phase = zeros(s.hh,ndft);
%out.Ht.amp = zeros(s.hh,ndft);
%out.Ht.phase = zeros(s.hh,ndft);
out.Hp.amp = zeros(s.hh,ndft);
out.Hp.phase = zeros(s.hh,ndft);

for m = 1:ndft,
	  Ertemp = dft.Er(2*m-1,:) + 1i*dft.Er(2*m,:);
out.Er.amp(:,m) = abs(Ertemp);
out.Er.phase(:,m) = unwrap(angle(Ertemp));
    
%Ettemp = dft.Et(2*m-1,:) + 1i*dft.Et(2*m,:);
%out.Et.amp(:,m) = abs(Ettemp);
%out.Et.phase(:,m) = unwrap(angle(Ettemp));
    
%Eptemp = dft.Ep(2*m-1,:) + 1i*dft.Ep(2*m,:);
%out.Ep.amp(:,m) = abs(Eptemp);
%out.Ep.phase(:,m) = unwrap(angle(Eptemp));
    
%Hrtemp = dft.Hr(2*m-1,:) + 1i*dft.Hr(2*m,:);
%out.Hr.amp(:,m) = abs(Hrtemp);
%out.Hr.phase(:,m) = unwrap(angle(Hrtemp));
    
%Httemp = dft.Ht(2*m-1,:) + 1i*dft.Ht(2*m,:);
%out.Ht.amp(:,m) = abs(Httemp);
%out.Ht.phase(:,m) = unwrap(angle(Httemp));
    
Hptemp = dft.Hp(2*m-1,:) + 1i*dft.Hp(2*m,:);
out.Hp.amp(:,m) = abs(Hptemp);
out.Hp.phase(:,m) = unwrap(angle(Hptemp));
    
end

out.dist = s.th*s.RE/1000;
out.dftfreqs = dftfreqs;

%% plotting

if (doplot),
    
  freqs = [5 10 15 20 25 30];
    
h1 = figure(1);
set(h1,'position',[500 200 1500 1200]);
for m = 1:ndft,
	  ax(m,1) = subplot(ndft,2,2*m-1);
ax(m,2) = subplot(ndft,2,2*m);
        
phasetoplot = out.Er.phase(1:1600,m) + (2*pi/3e8 * dftfreqs(m)*out.dist*1000);
plot(ax(m,1),out.dist,log10(out.Er.amp(1:1600,m)));
plot(ax(m,2),out.dist,phasetoplot);
    end
    
end
