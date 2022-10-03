%  Meyermeasurements_analyze.m
% 
% Peter Svensson 11 April 2022 & 20 April 2022

% Ivar's notes: 

indatafolder = 'Meyer4XPmeasurements_data/';
filenamestart = 'Meyer';
filenameend = '_S01_R01.etx';

filenumbers = [0:120];
nfiles = length(filenumbers);
numberofsampelstoread = 2000;

reloadfiles = 1;

nfft = 8192;
%ivfft = 320:1300;
ivfft = 320 + [0:479];  % 480 samples = 10 ms
ivfft = 320 + [0:239];  % 240 samples = 5 ms

plotalsosim = 0;

ivplot = 301:500;
shiftamplp = 0.2;
shiftamp = 0.5;

nfirlp = 101;
freqlpfir = 5000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if reloadfiles == 1
    disp(' ')
    for ii = 1:nfiles
        II = int2str(filenumbers(ii));
        if round(ii/10)*10 == ii
            disp(['   Loading file: ',II])
        end
        if filenumbers(ii) < 10
            II = ['00',II];
        elseif filenumbers(ii) < 100
            II = ['0',II];
        end
        filename = [indatafolder,filenamestart,II,filenameend];
        [ir,fs] = readetxfile(filename,numberofsampelstoread);
    
        if ii == 1
            allirs = zeros(length(ir),nfiles);
            allirslp = zeros(length(ir),nfiles);
            Blp = fir1(nfirlp,freqlpfir/(fs/2));
        end
        allirs(:,ii) = -ir;
        irf = contwo(-ir,Blp);
        irf = irf(floor(nfirlp/2):end);
        irf = irf(1:length(ir));
        allirslp(:,ii) = irf;

    %     figure(1)
    %     clf(1)
    %     plot(ir)
    %     grid
    %     pause
    
    end
end

shiftvec = [0:nfiles-1];
figure(1)
clf(1)
h = plot(allirslp(ivplot,:) + shiftamplp*shiftvec);
%h = plot(allirslp(ivplot,:));
seth
grid
g = title('All IRs, low-passfiltered at 5kHz, Meyer 4XP');
sg14
g = xlabel('Sample number   [-]');
sg14

figure(5)
clf(5)
h = plot(allirs(ivplot,:) + shiftamp*shiftvec);
seth
grid
g = title('All IRs, Meyer 4XP');
sg14
g = xlabel('Sample number   [-]');
sg14


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alltfs = fft(allirs(ivfft,:),nfft);

fvec = fs/nfft*[0:nfft/2-1];

figure(2)
clf(2)
h = semilogx(fvec,20*log10(abs(alltfs(1:nfft/2,:))));
seth
grid
g = xlabel('Frequency   [Hz]');
sg14
g = ylabel('FR magnitude   [dB]');
sg14
axis([80 20000 -10 50])
g = title('All frequency responses, Meyer 4XP');
sg14

intensity = abs(alltfs(1:nfft/2,1:end-1)).^2;
intensity_onaxis = intensity(:,1);
intensity_oneoff = intensity(:,2);
intensity_30 = intensity(:,1+10);
intensity_60 = intensity(:,1+20);
intensity_90 = intensity(:,1+30);
intensity_120 = intensity(:,1+40);
intensity_150 = intensity(:,1+50);
intensity_180 = intensity(:,1+60);
intensity_210 = intensity(:,1+70);
intensity_240 = intensity(:,1+80);
intensity_270 = intensity(:,1+90);
intensity_300 = intensity(:,1+100);
intensity_330 = intensity(:,1+110);
intensity_avg = mean(intensity.').';

figure(3)
clf(3)
h = semilogx(fvec,20*log10(abs(alltfs(1:nfft/2,[1 end]))),fvec,10*log10(intensity_avg));
%h = semilogx(fvec,20*log10(abs(alltfs(1:nfft/2,[1])))-20*log10(abs(alltfs(1:nfft/2,[end]))));
seth
grid
for ii = 1:length(h)
    set(h(ii),'LineWidth',2)
end
g = xlabel('Frequency   [Hz]');
sg14
g = ylabel('FR magnitude   [dB]');
sg14
g = title('Frequency response, Meyer 4XP');
sg14
g = legend('0 deg.','360 deg.','Average (2D)');
sg14
set(g,'Location','SouthEast');
xlim([80 20000])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DF_2D = intensity_onaxis./intensity_avg;
DI_2D = 10*log10(DF_2D);

DF_2D_30 = intensity_30./intensity_avg;
DI_2D_30 = 10*log10(DF_2D_30);

DF_2D_60 = intensity_60./intensity_avg;
DI_2D_60 = 10*log10(DF_2D_60);

DF_2D_90 = intensity_90./intensity_avg;
DI_2D_90 = 10*log10(DF_2D_90);

DF_2D_120 = intensity_120./intensity_avg;
DI_2D_120 = 10*log10(DF_2D_120);

DF_2D_150 = intensity_150./intensity_avg;
DI_2D_150 = 10*log10(DF_2D_150);

DF_2D_180 = intensity_180./intensity_avg;
DI_2D_180 = 10*log10(DF_2D_180);

DF_2D_210 = intensity_210./intensity_avg;
DI_2D_210 = 10*log10(DF_2D_210);

DF_2D_240 = intensity_240./intensity_avg;
DI_2D_240 = 10*log10(DF_2D_240);

DF_2D_270 = intensity_270./intensity_avg;
DI_2D_270 = 10*log10(DF_2D_270);

DF_2D_300 = intensity_300./intensity_avg;
DI_2D_300 = 10*log10(DF_2D_300);

DF_2D_330 = intensity_330./intensity_avg;
DI_2D_330 = 10*log10(DF_2D_330);

DF_2D_oneoff = intensity_oneoff./intensity_avg;
DI_2D_oneoff = 10*log10(DF_2D_oneoff);

figure(4)
if plotalsosim == 1
%    h = semilogx(fvec_tf,DI_sim,fvec,DI_2D,fvec,DI_2D_oneoff,fvec_tf,DI_2D_sim);
    h = semilogx(fvec,DI_2D+15,   fvec_tf,DI_2D_sim+15,'--',...
                 fvec,DI_2D_30+10,fvec,DI_2D_330+10,fvec_tf,DI_2D_sim_30+10,'--',...
                 fvec,DI_2D_60+5,fvec,DI_2D_300+5,fvec_tf,DI_2D_sim_60+5,'--',...
                 fvec,DI_2D_90,fvec,DI_2D_270,fvec_tf,DI_2D_sim_90,'--',...     
                 fvec,DI_2D_120-5,fvec,DI_2D_240-5,fvec_tf,DI_2D_sim_120-5,'--',...
                 fvec,DI_2D_150-10,fvec,DI_2D_210-10,fvec_tf,DI_2D_sim_150-10,'--',...
                 fvec,DI_2D_180-15,fvec_tf,DI_2D_sim_180-15,'--');
else
    h = semilogx(fvec,DI_2D);
end
%seth
for ii = 1:length(h)
    if ismember(ii,[1 3 6 9 12 15 18])
        linecol = get(h(ii),'Color');
    else
        set(h(ii),'Color',linecol);
    end
    set(h(ii),'LineWidth',1)
end
grid
g = xlabel('Frequency   [Hz]');
sg14
g = ylabel('DI   [dB]');
sg14
g = title('Hor. plane DI for Meyer 4XP, at 2.35m distance');
sg14
axis([80 20000 -35 30])
if plotalsosim == 1
    g = legend('Meas., on-axis +15dB','Sim., on-axis +15dB',...
        'Meas., 30 deg. +10dB','Meas., 330 deg. +10dB','Sim., 30 deg. +10dB',...
        'Meas., 60 deg. +5dB','Meas., 300 deg. +5dB','Sim., 60 deg. +5dB',...
        'Meas., 90 deg.','Meas., 270 deg.','Sim., 90 deg.',...
        'Meas., 120 deg. -5dB','Meas., 240 deg. -5dB','Sim., 120 deg. -5dB',...
        'Meas., 150 deg. -10dB','Meas., 210 deg. -10dB','Sim., 150 deg. -10dB',...
        'Meas., 180 deg. -15dB','Sim., 180 deg. -15dB');
else
    g = legend('Measured, 2D DI');
end
set(g,'Location','Eastoutside','FontSize',12)


[LIA,LOCB] = ismember(fvec,fvec_tf);
iv_fvec = find(LIA);

figure(5)
clf(5)
h = semilogx(fvec(iv_fvec),DI_2D(iv_fvec)-DI_2D_sim,'-',...
             fvec(iv_fvec),DI_2D_30(iv_fvec)-DI_2D_sim_30,'-',...
             fvec(iv_fvec),DI_2D_330(iv_fvec)-DI_2D_sim_30,'--',...
             fvec(iv_fvec),DI_2D_60(iv_fvec)-DI_2D_sim_60,'-',...
             fvec(iv_fvec),DI_2D_300(iv_fvec)-DI_2D_sim_60,'--',...
             fvec(iv_fvec),DI_2D_90(iv_fvec)-DI_2D_sim_90,'-',...
             fvec(iv_fvec),DI_2D_270(iv_fvec)-DI_2D_sim_90,'--',...
             fvec(iv_fvec),DI_2D_120(iv_fvec)-DI_2D_sim_120,'-',...
             fvec(iv_fvec),DI_2D_240(iv_fvec)-DI_2D_sim_120,'--',...
             fvec(iv_fvec),DI_2D_150(iv_fvec)-DI_2D_sim_150,'-',...
             fvec(iv_fvec),DI_2D_210(iv_fvec)-DI_2D_sim_150,'--',...
             fvec(iv_fvec),DI_2D_180(iv_fvec)-DI_2D_sim_180,'-');
seth
grid
for ii = 1:length(h)
    if ismember(ii,[1 2 4 6 8 10 12])
        linecol = get(h(ii),'Color');
    else
        set(h(ii),'Color',linecol);
    end
    set(h(ii),'LineWidth',1)
end
g = xlabel('Frequency   [Hz]');
sg14
g = ylabel('DI-difference   [dB]');
sg14
g = title('Hor. plane DI, measured - simulated for Meyer 4XP, at 2.35m distance');
sg14
axis([100 1500 -1 1])
g = legend('On-axis','30 deg.','330 deg.',...
                     '60 deg.','300 deg.',...
                     '90 deg.','270 deg.',...
                     '120 deg.','240 deg.',...
                     '150 deg.','210 deg.',...
                     '180 deg.');
set(g,'Location','Eastoutside','FontSize',12)


ivhist = find(fvec(iv_fvec) > 100 & fvec(iv_fvec) < 1500);
diffvals = DI_2D(iv_fvec(ivhist))-DI_2D_sim(ivhist);
diffvals = [diffvals;DI_2D_30(iv_fvec(ivhist))-DI_2D_sim_30(ivhist)];
diffvals = [diffvals;DI_2D_60(iv_fvec(ivhist))-DI_2D_sim_60(ivhist)];
diffvals = [diffvals;DI_2D_90(iv_fvec(ivhist))-DI_2D_sim_90(ivhist)];
diffvals = [diffvals;DI_2D_120(iv_fvec(ivhist))-DI_2D_sim_120(ivhist)];
diffvals = [diffvals;DI_2D_150(iv_fvec(ivhist))-DI_2D_sim_150(ivhist)];
diffvals = [diffvals;DI_2D_180(iv_fvec(ivhist))-DI_2D_sim_180(ivhist)];
diffvals = [diffvals;DI_2D_210(iv_fvec(ivhist))-DI_2D_sim_150(ivhist)];
diffvals = [diffvals;DI_2D_240(iv_fvec(ivhist))-DI_2D_sim_120(ivhist)];
diffvals = [diffvals;DI_2D_270(iv_fvec(ivhist))-DI_2D_sim_90(ivhist)];
diffvals = [diffvals;DI_2D_300(iv_fvec(ivhist))-DI_2D_sim_60(ivhist)];
diffvals = [diffvals;DI_2D_330(iv_fvec(ivhist))-DI_2D_sim_30(ivhist)];
diffvals_LF = diffvals;

xvals = diffvals_LF(44*1+1:44*2);
yvals = diffvals_LF(44*11+1:44*12);
xvals = [xvals;diffvals_LF(44*2+1:44*3)];
yvals = [yvals;diffvals_LF(44*10+1:44*11)];
xvals = [xvals;diffvals_LF(44*3+1:44*4)];
yvals = [yvals;diffvals_LF(44*9+1:44*10)];
xvals = [xvals;diffvals_LF(44*4+1:44*5)];
yvals = [yvals;diffvals_LF(44*8+1:44*9)];
xvals = [xvals;diffvals_LF(44*5+1:44*6)];
yvals = [yvals;diffvals_LF(44*7+1:44*8)];

[rho,pval] = corr([xvals yvals]);


ivhist = find(fvec(iv_fvec) > 2000);
diffvals = DI_2D(iv_fvec(ivhist))-DI_2D_sim(ivhist);
diffvals = [diffvals;DI_2D_30(iv_fvec(ivhist))-DI_2D_sim_30(ivhist)];
diffvals = [diffvals;DI_2D_60(iv_fvec(ivhist))-DI_2D_sim_60(ivhist)];
diffvals = [diffvals;DI_2D_90(iv_fvec(ivhist))-DI_2D_sim_90(ivhist)];
diffvals = [diffvals;DI_2D_120(iv_fvec(ivhist))-DI_2D_sim_120(ivhist)];
diffvals = [diffvals;DI_2D_150(iv_fvec(ivhist))-DI_2D_sim_150(ivhist)];
diffvals = [diffvals;DI_2D_180(iv_fvec(ivhist))-DI_2D_sim_180(ivhist)];
diffvals = [diffvals;DI_2D_210(iv_fvec(ivhist))-DI_2D_sim_150(ivhist)];
diffvals = [diffvals;DI_2D_240(iv_fvec(ivhist))-DI_2D_sim_120(ivhist)];
diffvals = [diffvals;DI_2D_270(iv_fvec(ivhist))-DI_2D_sim_90(ivhist)];
diffvals = [diffvals;DI_2D_300(iv_fvec(ivhist))-DI_2D_sim_60(ivhist)];
diffvals = [diffvals;DI_2D_330(iv_fvec(ivhist))-DI_2D_sim_30(ivhist)];
diffvals_HF = diffvals;

figure(6)
hist(diffvals_LF,20);
g = xlabel('Measurement - simulation   [dB]');
sg14
g = ylabel('No. of values   [-]');
sg14
grid

quantvals_LF = quantile(diffvals_LF,[0.025 0.05 0.5 0.95 0.975]);










return

figure(8)
clf(8)
df = fs/nfft;
freqtocheck = 200:10:1e4;
angshift = zeros(size(freqtocheck));
measanglevec = [0:3:357];
iv1 = find(measanglevec>0 & measanglevec<180);
iv2 = find(measanglevec>180 & measanglevec<360);

for ii = 1:length(freqtocheck)
    ivfindfreq = round(freqtocheck(ii)/df)+1;
    F1 = abs(alltfs(ivfindfreq,iv1));
    F2 = abs(alltfs(ivfindfreq,iv2(end:-1:1)));
    F1 = F1(:);
    F2 = F2(:);
    maxval = max(max([F1 F2]));
    minval = min(min([F1 F2]));
    halfval = minval + 0.6*(maxval-minval);
    iv1half = find(F1<halfval);
    iv1half = iv1half(1)
    iv2half = find(F2<halfval);
    iv2half = iv2half(1);
    ang1 = measanglevec(iv1(iv1half));
    ang2 = measanglevec(iv1(iv2half));
    angshift(ii) = ang1 - ang2;

%     h = plot([1:length(F1)],[F1(:) F2(:)],'-o',iv1half,F1(iv1half),'b*',iv2half,F2(iv2half),'r*');
%     seth
%     grid
%     g = xlabel('Radiation angle no.  [-]');
%     %sg14
%     g = ylabel('Frequency response mag.   [-]');
%     %sg14
%     legend([num2str(ang1),' deg.'],[num2str(ang2),' deg.'])
%     title(['Freq.: ',num2str(fvec(ivfindfreq)),' Hz'])
%     pause
end

figure(9)
clf(9)
h = semilogx(freqtocheck,angshift,'-o');
seth
grid
g = xlabel('Frequency   [Hz]');
sg14
g = ylabel('Angle shift   [deg.]');
sg14



