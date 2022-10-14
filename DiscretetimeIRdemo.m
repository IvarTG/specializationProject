% Demo of discretization error for impulse responses

r1 = 0.5;
r2 = 0.55;

% samplemethod = 1 -> round at low fs
% samplemethod = 2 -> round at high fs and decimate
% samplemethod = 3 -> round at low fs but split pulse between two samples

samplemethod = 1;

% Generate discrete-time IRs
% Exact arrival times
cair = 343;
t1 = r1/cair;
t2 = r2/cair;

% Approach 1: round to nearest discrete time-sample
fs = 48000;

if samplemethod == 2
    noversamp = 8;

    fs_oversamp = fs*noversamp;
    dt = 1/fs_oversamp;
elseif samplemethod == 1 
    dt = 1/fs;
elseif samplemethod == 3 
    dt = 1/fs;
end
n1_exact = t1/dt + 1;
n2_exact = t2/dt + 1;

n1_rounded = round(t1/dt) + 1;
n2_rounded = round(t2/dt) + 1;

irmaxlength = max([n1_rounded n2_rounded]) + 1;
ir1 = zeros(irmaxlength+500,1);
ir2 = zeros(irmaxlength+500,1);

if samplemethod < 3
    ir1(n1_rounded) = 1/r1; %IR(t1) ~= IR(n1_rounded/fs) = dirac(t - r1/c)/r1 = dirac((fs*r1)/(c*fs) - r1/c)/r1 = dirac(0)/r1 = 1/r1??
    ir2(n2_rounded) = 1/r2;
else
    n1_frac = n1_exact - floor(n1_exact);
    n2_frac = n2_exact - floor(n2_exact);

    ir1(floor(n1_exact)) = 1/r1*(1-n1_frac);
    ir1(floor(n1_exact)+1) = 1/r1*n1_frac;
    ir2(floor(n2_exact)) = 1/r2*(1-n2_frac);
    ir2(floor(n2_exact)+1) = 1/r2*n2_frac;
end
% figure(3)
% plot(ir1)
% pause

if samplemethod == 2
    ir1 = noversamp*decimate(ir1,noversamp);
    ir2 = noversamp*decimate(ir2,noversamp);
end

% figure(3)
% plot([ir1 ir2])
% pause

nfft = 8192;
F1 = fft(ir1,nfft);
F2 = fft(ir2,nfft);

fvec_fft = fs/nfft*[0:nfft/2-1].';

ptot_fromIR = F1 + F2;
ptot_fromIR = ptot_fromIR(1:nfft/2);

% Truth: FD results

cair = 343;
kvec = 2*pi*fvec_fft/cair;

% Correct answer, in freq. domain
ptot = exp(-1i*kvec*r1)/r1 + exp(-1i*kvec*r2)/r2;

if samplemethod == 2
    % Correct scale for IR results
    %ptot_fromIR = ptot_fromIR*(1/r1 + 1/r2)/(sum(ir1)+sum(ir2));
end

figure(1)
clf(1)
h = plot(fvec_fft,abs(ptot),'o',fvec_fft,abs(ptot_fromIR));
seth
grid
legend('FD results (truth)','IR-to-FR results');

figure(3)
clf(3)
h = plot(fvec_fft,unwrap(angle(ptot)),'o',fvec_fft,unwrap(angle(ptot_fromIR)));
seth
grid
legend('FD results (truth)','IR-to-FR results');
g = ylabel('Phase   [rad].');
sg14

resultdiff = ptot_fromIR(:) - ptot(:);
figure(2)
h = plot(fvec_fft,abs(resultdiff)./abs(ptot),'-');
seth
grid
g = ylabel('Rel. error   [-]');
sg14












