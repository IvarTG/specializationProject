%  Meyermeasurements_analyze.m
% 
% Peter Svensson 11 April 2022 & 20 April 2022

indatafolder = 'Meyer4XPmeasurements_data/';
filenamestart = 'Meyer';
filenameend = '_S01_R01.etx';

filenumbers = [0:120];
nfiles = length(filenumbers);
numberofsampelstoread = 2000;

reloadfiles = 0;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Discrete IR Code
%clf
r = 2.35; %m
xc = 0;   %center x
yc = 0;   %center y

n_meas = 24;

theta = linspace(0,2*pi-((3*pi)/180),n_meas);
x = r*cos(theta) + xc;
y = r*sin(theta) + yc;

mic_xyz = cat(1, x, y, zeros(1, n_meas));

ES1 = [0,0,0];
ES2 = [0.3,0,0];
ES3 = [-0.3,0,0];

ES_xyz = cat(1, ES1, ES2, ES3);
n_source = 3;

figure(2)
scatter(x,y)
hold on
plot(ES1(1),ES1(2), 'x', 'LineWidth', 2, 'MarkerSize', 9);
plot(ES2(1),ES2(2), 'x', 'LineWidth', 2, 'MarkerSize', 9);
plot(ES3(1),ES3(2), 'x', 'LineWidth', 2, 'MarkerSize', 9);
grid
axis equal
hold off

r_ES = EDcalcdist(ES_xyz, mic_xyz.');
r_ES1 = r_ES(1,:);
r_ES2 = r_ES(2,:);
r_ES3 = r_ES(3,:);

% Discretization error for impulse responses

% samplemethod = 1 -> round at low fs
% samplemethod = 2 -> round at high fs and decimate
% samplemethod = 3 -> round at low fs but split pulse between two samples

samplemethod = 2;

% Generate discrete-time IRs
% Exact arrival times
cair = 343;

% Approach 1: round to nearest discrete time-sample
fs = 48000;

if samplemethod == 2
    noversamp = 8;
    fs_oversamp = fs*noversamp;
    dt = 1/fs_oversamp;
elseif samplemethod == 1 % Approach 1: round to nearest discrete time-sample
    dt = 1/fs;
elseif samplemethod == 3 
    dt = 1/fs;
end

n_exact = zeros(3, n_meas);
n_rounded = zeros(3, n_meas);

for i = 1:n_source
    for j = 1:n_meas
        t_ij = r_ES(i,j)/cair;
        n_exact(i,j) = t_ij/dt + 1;
        n_rounded(i,j) = round(t_ij/dt) + 1;
    end
end

irmaxlength = -inf;
for i = 1:n_source
    if max(n_rounded(i,:)) > irmaxlength
        irmaxlength = max(n_rounded(i,:));
    end
end
irmaxlength = irmaxlength+1;

%irmaxlength = max([n_rounded(1,:) n_rounded(2,:) n_rounded(3,:)]) + 1;
ir = zeros(irmaxlength+5000, n_source, n_meas);
nfft = 8192;
F = zeros(nfft,n_source, n_meas);

if samplemethod == 1
    for i = 1:n_source
        for j = 1:n_meas
            ir(n_rounded(i,j),i,j) = 1/r_ES(i,j);
            F(:,i,j) = fft(ir(:,i,j),nfft);
        end
    end
elseif samplemethod == 2
    c = zeros(ceil((irmaxlength+5000)/noversamp), n_source, n_meas);
    for i = 1:n_source
        for j = 1:n_meas
            ir(n_rounded(i,j),i,j) = 1/r_ES(i,j);
            w = ir(:,i,j);
            c(:,i,j) = noversamp*decimate(w,noversamp);
            F(:,i,j) = fft(c(:,i,j),nfft);
        end
    end
    ir = c;
elseif samplemethod == 3
    for i = 1:n_source
        for j = 1:n_meas
            frac = n_exact(i,j) - floor(n_exact(i,j));
            ir(floor(n_exact(i,j)),i,j) = 1/r_ES(i,j)*(1-frac);
            ir(floor(n_exact(i,j))+1,i,j) = 1/r_ES(i,j)*frac;
            
            F(:,i,j) = fft(ir(:,i,j),nfft);
        end
    end
end

% if samplemethod == 2
%     ir = noversamp*decimate(ir,noversamp);
%     ir = reshape(ir,[],3,n_meas);
% %                 c = noversamp*decimate(squeeze(ir(i,j,:)),noversamp);
% %                 c = reshape(c,1,1,[]);
% %                 ir(i,j,:) = c;
% end

% figure(3)
% plot(ir(1,1))
% pause

fvec_fft = fs/nfft*[0:nfft/2-1].';

ptot_fromIR = sum(F,[2 3]);
ptot_fromIR = ptot_fromIR(1:nfft/2);
ptot_fromIR = squeeze(ptot_fromIR);

% Truth: FD results

cair = 343;
kvec = 2*pi*fvec_fft/cair;

% Correct answer, in freq. domain
ptot = 0;
for i = 1:n_source
    for j = 1:n_meas
        ptot = ptot + exp(-1i*kvec*r_ES(i,j))/r_ES(i,j);
    end     
end
if samplemethod == 2
    % Correct scale for IR results
    %ptot_fromIR = ptot_fromIR*(1/r1 + 1/r2)/(sum(ir1)+sum(ir2));
    %ptot_fromIR = ptot_fromIR*(sum(1./r1))/(sum(ir, 'all'));
end

