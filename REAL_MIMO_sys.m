clf
close all
%clear

indatafolder = 'Meyer4XPmeasurements_data/';
filenamestart = 'Meyer';
filenameend = '_S01_R01.etx';

filenumbers = [0:120];
nfiles = length(filenumbers);
numberofsampelstoread = 2000;

reloadfiles =1;

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

if reloadfiles == 0
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
    end
end

%Discrete IR Code
%clf
r = 2.35; %m
xc = 0;   %center x
yc = 0;   %center y

n_meas = 24;

theta = linspace(0,2*pi-((3*pi)/180),n_meas); %if n_meas is 120
%theta = linspace(0,2*pi,n_meas);                %if n_meas is 121
x = r*cos(theta) + xc;
y = r*sin(theta) + yc;

mic_xyz = cat(1, x, y, zeros(1, n_meas));

ES1 = [0.0,0,0];
ES2 = [0.03,0,0];
ES3 = [0.25,0,0];
% ES1 = [0.05,-0.05,0];
% ES2 = [0.07,0.1,0];
% ES3 = [0.07,-0.1,0];
ES4 = [0.07,0,0];
ES5 = [0.12,0,0];


ES_xyz = cat(1, ES1, ES2, ES3, ES4, ES5);
n_source = 5;

figure(2)
scatter(x,y)
hold on
plot(ES1(1),ES1(2), 'x', 'LineWidth', 2, 'MarkerSize', 9);
plot(ES2(1),ES2(2), 'x', 'LineWidth', 2, 'MarkerSize', 9);
plot(ES3(1),ES3(2), 'x', 'LineWidth', 2, 'MarkerSize', 9);
plot(ES4(1),ES4(2), 'x', 'LineWidth', 2, 'MarkerSize', 9);
plot(ES5(1),ES5(2), 'x', 'LineWidth', 2, 'MarkerSize', 9);
grid
axis equal
hold off

r_ES = EDcalcdist(ES_xyz, mic_xyz.');
r_ES1 = r_ES(1,:);
r_ES2 = r_ES(2,:);
r_ES3 = r_ES(3,:);
r_ES = r_ES(1:5,:);

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

n_exact = zeros(n_source, n_meas);
n_rounded = zeros(n_source, n_meas);

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

ir = zeros(irmaxlength+2500, n_source, n_meas);
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
    c = zeros(ceil((irmaxlength+2500)/noversamp), n_source, n_meas);
    for i = 1:n_source
        for j = 1:n_meas
            frac = n_exact(i,j) - floor(n_exact(i,j));
            ir(floor(n_exact(i,j)),i,j) = 1/r_ES(i,j)*(1-frac);
            ir(floor(n_exact(i,j))+1,i,j) = 1/r_ES(i,j)*frac;
            w = ir(:,i,j);
            c(:,i,j) = noversamp*decimate(w,noversamp,1000,'fir');
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

xinput = ir;
n_steps = length(xinput);
%h_length = 2e2;
h_length = length(ir);
h_startestimate = zeros(h_length,n_source);
y_step = 120/n_meas;
% youtput = cat(1, allirs(:,1:y_step:120), zeros(3000,n_meas));
youtput = allirslp(:,1:y_step:120);
w1 = window(@hamming,135);
% test2 = zeros(length(youtput),n_meas);
% test2(290:424,:) = youtput(290:424,:).*w1;
% youtput = test2;
% convfactor = 0.00009;
convfactor = 0.0001;
%convfactor = 0.00234;

[q_est,errorhistory] = MC_LMS_Ncoeffs(xinput,youtput,n_source,n_meas,h_length,convfactor,n_steps);
totalerrorhistory = errorhistory(end,:);

i = 22;

while i > 0
    [q_est,errorhistory] = MC_LMS_Ncoeffs(xinput,youtput,n_source,n_meas,h_length,convfactor,n_steps,q_est);
    totalerrorhistory = [totalerrorhistory; errorhistory(end,:)];
    i = i-1;
end
figure(7)
semilogy(abs(totalerrorhistory))
title('errorhist')
grid


y_est_2_len = contwo(q_est, ir);
y_est_2 = zeros(length(y_est_2_len), n_meas);
for i = 1:n_meas
    for j = 1:n_source
        y_est_2(:,i) = y_est_2(:,i) + contwo(q_est(:,j), ir(:,j,i));
    end
end

y_est_1 = zeros(length(youtput), n_meas);
for ii =1:n_steps
    for kk = 1:n_meas
        for jj = 1:n_source
        
            flippedx = ir(ii:-1:max([1 ii-h_length+1]),jj,kk);
            lengthflippedx = length(flippedx);
            if lengthflippedx < h_length
              flippedx = [flippedx;zeros(h_length-lengthflippedx,1)];
            end
            y_est_1(ii,kk) = y_est_1(ii,kk) + sum(flippedx.*q_est(:,jj));
            
        end
    end
end

% for ii =1:length(q_est)
%     for jj = 1:n_source
%         for kk = 1:n_meas
%             flippedx = q_est(ii:-1:max([1 ii-h_length+1]),jj);
%             lengthflippedx = length(flippedx);
%             if lengthflippedx < h_length
%               flippedx = [flippedx;zeros(h_length-lengthflippedx,1)];
%             end
%             y_est_1(ii,kk) = y_est_1(ii,kk) + sum(flippedx.*ir(:,jj,kk));
%             
%         end
%     end
% end

nfft = 4096;
fvec = fs/nfft*[0:nfft/2-1];
froutput = fft(youtput(:,1),nfft); 
froutput = froutput(1:nfft/2);
frest = fft(y_est_1(:,1),nfft); 
frest = frest(1:nfft/2);

figure(9)
% plot(lowpass(y_est_1(:,1),7000,48000))
plot(y_est_1(:,1),'-k')
xlim([200 800])
hold on
plot(y_est_2(:,1),'*')
xlim([200 800])
title('y_{est} vs y_{real}')
plot(youtput(:,1),'--')
legend('y_{est1}','y_{est1}','y_{real}')

figure(10)
semilogx(fvec,20*log10(abs(froutput)),fvec,20*log10(abs(frest)),'o' )
legend('froutput','frest')
grid

figure(11)
plot(y_est_1)
title('y_{est1}')
figure(19)
plot(y_est_2)
title('y_{est2}')
xlim([200 600])

figure(12)
plot(y_est_1+0.2*[0:n_meas-1])
xlim([200 600])
figure(13)
plot(youtput + 0.2*[0:n_meas-1])
xlim([200 600])

figure(14)
subplot(4,1,1)
plot(q_est(:,1))
title('q_{est1}')
subplot(4,1,2)
plot(q_est(:,2))
title('q_{est2}')
% subplot(4,1,3)
% plot(q_est(:,3))
% title('q_{est3}')
% subplot(4,1,4)
% plot(q_est(:,4))
% title('q_{est4}')