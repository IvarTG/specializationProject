clf
close all
clear

%Discrete IR Code
%clf
r = 2.35; %m
xc = 0;   %center x
yc = 0;   %center y

n_meas = 24;

%theta = linspace(0,2*pi-((3*pi)/180),n_meas); %if n_meas is 120
theta = linspace(0,2*pi,n_meas);                %if n_meas is 121
x = r*cos(theta) + xc;
y = r*sin(theta) + yc;

mic_xyz = cat(1, x, y, zeros(1, n_meas));

ES1 = [0,0,0];
ES2 = [-0.32,0,0];
ES3 = [-0.3,0,0];

ES_xyz = cat(1, ES1, ES2, ES3);
n_source = 1;

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
r_ES = r_ES2;

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

%irmaxlength = max([n_rounded(1,:) n_rounded(2,:) n_rounded(3,:)]) + 1;
ir = zeros(irmaxlength+500, n_source, n_meas);
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
    c = zeros(ceil((irmaxlength+500)/noversamp), n_source, n_meas);
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

% Simulated Input q_sim
q_sim = zeros(length(ir),1);
q_sim(69) = 1;
%q_sim = randn(1e2,1);
%q_sim = 1;

%Simulated Output y_sim
n_steps = length(q_sim);
y_sim_1 = zeros(n_steps,n_meas);
h_length = length(ir);

for ii =1:length(q_sim)
    for jj = 1:n_source
        for kk = 1:n_meas
            flippedx = q_sim(ii:-1:max([1 ii-h_length+1]),jj);
            lengthflippedx = length(flippedx);
            if lengthflippedx < h_length
              flippedx = [flippedx;zeros(h_length-lengthflippedx,1)];
            end
            y_sim_1(ii,kk) = y_sim_1(ii,kk) + sum(flippedx.*ir(:,jj,kk));
            
        end
    end
end


y_sim_2_len = contwo(q_sim, ir);
y_sim_2 = zeros(length(y_sim_2_len)+500, n_meas);

for i = 1:n_meas
    for j = 1:n_source
        y_sim_2(:,i) = y_sim_2(:,i) + sum(contwo(q_sim(:,j), ir(:,j,i)));
    end
end

xinput = ir;
n_steps = length(xinput);
h_length = 2e2;
% h_length = length(ir);
h_startestimate = zeros(h_length,n_source);
youtput = y_sim_1;
convfactor = 0.01;

[q_est,errorhistory] = MC_LMS_Ncoeffs(xinput,youtput,n_source,n_meas,h_length,convfactor,n_steps);
totalerrorhistory = errorhistory(end,:);

i = 10;

while i > 0
    [q_est,errorhistory] = MC_LMS_Ncoeffs(xinput,youtput,n_source,n_meas,h_length,convfactor,n_steps,q_est);
    totalerrorhistory = [totalerrorhistory; errorhistory(end,:)];
    i = i-1;
end
figure(7)
semilogy(abs(totalerrorhistory))
title('errorhist')
grid

 figure(8)
 plot(q_sim,'k');
 title('q_{sim} vs. q_{est}')
 hold on
 grid
 plot(q_est,'--')
 legend('q_{sim}', 'q_{est}')

y_est_2_len = contwo(q_est, ir);
y_est_2 = zeros(length(y_est_2_len), n_meas);
for i = 1:n_meas
    for j = 1:n_source
        y_est_2(:,i) = y_est_2(:,i) + contwo(q_est(:,j), ir(:,j,i));
    end
end

y_est_1 = zeros(length(youtput), n_meas);
for ii =1:n_steps
    for jj = 1:n_source
        for kk = 1:n_meas
            flippedx = ir(ii:-1:max([1 ii-h_length+1]),jj,kk);
            lengthflippedx = length(flippedx);
            if lengthflippedx < h_length
              flippedx = [flippedx;zeros(h_length-lengthflippedx,1)];
            end
            y_est_1(ii,kk) = y_est_1(ii,kk) + sum(flippedx.*q_est(:,jj));
            
        end
    end
end
nfft = 4096;
fvec = fs/nfft*[0:nfft/2-1];
froutput = fft(youtput(:,1),nfft); 
froutput = froutput(1:nfft/2);
frest = fft(y_est_1(:,1),nfft); 
frest = frest(1:nfft/2);

figure(9)
% plot(lowpass(y_est_1(:,1),7000,48000))
plot(y_est_1(:,20),'-k')
xlim([200 800])
hold on
title('y_{est} vs y_{real}')
plot(youtput(:,20),'--')
legend('y_{est}','y_{real}')

figure(10)
semilogx(fvec,20*log10(abs(froutput)),fvec,20*log10(abs(frest)),'o' )
legend('froutput','frest')
grid

figure(11)
plot(y_est_1)
title('y_{est1}')
figure(12)
plot(y_est_2)
title('y_{est2}')
%xlim([200 600])
figure(13)
plot(y_sim_1)
title('y_{sim1}')
%xlim([200 600])
figure(14)
plot(y_sim_2)
title('y_{sim2}')
%xlim([200 600])