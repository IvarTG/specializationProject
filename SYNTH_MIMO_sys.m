clf
close all
clear
%% Find 
r = 2.35; %m
xc = 0;   %center x
yc = 0;   %center y

n_meas = 24;

theta = linspace(0,2*pi-((3*pi)/180),n_meas); %if n_meas is 120
%theta = linspace(0,2*pi,n_meas);                %if n_meas is 121
x = r*cos(theta) + xc;
y = r*sin(theta) + yc;

mic_xyz = cat(1, x, y, zeros(1, n_meas));

OS1 = [0.0,0,0];
OS2 = [0.3,0,0];
OS3 = [0.25,0,0];

OS_xyz = cat(1, OS1, OS2, OS3);
n_source_OS = 3;

figure(2)
scatter(x,y)
hold on
plot(OS1(1),OS1(2), 'x', 'LineWidth', 2, 'MarkerSize', 9);
plot(OS2(1),OS2(2), 'x', 'LineWidth', 2, 'MarkerSize', 9);
plot(OS3(1),OS3(2), 'x', 'LineWidth', 2, 'MarkerSize', 9);
grid
axis equal
hold off

r_OS = EDcalcdist(OS_xyz, mic_xyz.');
r_OS1 = r_OS(1,:);
r_OS2 = r_OS(2,:);
r_OS3 = r_OS(3,:);
% r_ES_OS = r_ES_OS(1:2,:);

samplemethod = 2;

% Generate discrete-time IRs
% Exact arrival times
cair = 343;

% Approach 1: round to nearest discrete time-sample
fs = 48000;

ir_OS = find_IR(r_OS, n_meas,n_source_OS,samplemethod,cair,fs);

% Simulated Input q_sim
q_OS = zeros(length(ir_OS),n_source_OS);
q_OS(20,1) = 1;
q_OS(22,2) = 0.6;
q_OS(16,3) = 1.2;

%Simulated Output y_sim
n_steps = length(q_OS);
y_OS = zeros(n_steps,n_meas);
h_length = length(ir_OS);

for ii =1:length(q_OS)
    for jj = 1:n_source_OS
        for kk = 1:n_meas
            flippedx = ir_OS(ii:-1:max([1 ii-h_length+1]),jj,kk);
            lengthflippedx = length(flippedx);
            if lengthflippedx < h_length
              flippedx = [flippedx;zeros(h_length-lengthflippedx,1)];
            end
            y_OS(ii,kk) = y_OS(ii,kk) + sum(flippedx.*q_OS(:,jj));
            
        end
    end
end

%Find the IR signal of the estimatet sources q_ES

ES1 = [0.0,0,0];
ES2 = [0.3,0,0];
ES3 = [0.25,0,0];
ES4 = [0.4,0,0];
ES5 = [0.35,0,0];

%ES_xyz = cat(1, ES1, ES2, ES3);
ES_xyz = cat(1, ES1, ES2, ES3, ES4, ES5);
n_source_ES = 5;

r_ES = EDcalcdist(ES_xyz, mic_xyz.');

ir_ES = find_IR(r_ES, n_meas,n_source_ES,samplemethod,cair,fs);

xinput = ir_ES;
n_steps = length(xinput);
h_length = 2e2;
% h_length = length(ir);
h_startestimate = zeros(h_length,n_source_OS);
%youtput = y_OS;
youtput = cat(1, y_OS, zeros(length(xinput)-length(y_OS),n_meas));
%youtput = youtput + 0.001*randn(size(youtput));
convfactor = 0.00114;

[q_ES,errorhistory] = MC_LMS_Ncoeffs(xinput,youtput,n_source_ES,n_meas,h_length,convfactor,n_steps);
totalerrorhistory = errorhistory(end,:);

i = 60;

while i > 0
    [q_ES,errorhistory] = MC_LMS_Ncoeffs(xinput,youtput,n_source_ES,n_meas,h_length,convfactor,n_steps,q_ES);
    totalerrorhistory = [totalerrorhistory; errorhistory(end,:)];
    i = i-1;
end
figure(7)
semilogy(abs(totalerrorhistory))
title('errorhist')
grid

 figure(8)
 plot(q_OS(:,1),'k');
 title('q_{sim} vs. q_{est}')
 hold on
 grid
 plot(q_ES(:,1),'--')
 legend('q_{sim}', 'q_{est}')
 
 figure(18)
 plot(q_OS(:,2),'k');
 title('q_{sim} vs. q_{est}')
 hold on
 grid
 plot(q_ES(:,2),'--')
 legend('q_{sim}', 'q_{est}')
 
 figure(28)
 plot(q_OS(:,3),'k');
 title('q_{sim} vs. q_{est}')
 hold on
 grid
 plot(q_ES(:,3),'--')
 legend('q_{sim}', 'q_{est}')

y_est_2_len = contwo(q_ES, ir_ES);
y_est_2 = zeros(length(y_est_2_len), n_meas);
for i = 1:n_meas
    for j = 1:n_source_OS
        y_est_2(:,i) = y_est_2(:,i) + contwo(q_ES(:,j), ir_ES(:,j,i));
    end
end

y_est_1 = zeros(length(youtput), n_meas);
for ii =1:n_steps
    for kk = 1:n_meas
        for jj = 1:n_source_ES
        
            flippedx = ir_ES(ii:-1:max([1 ii-h_length+1]),jj,kk);
            lengthflippedx = length(flippedx);
            if lengthflippedx < h_length
              flippedx = [flippedx;zeros(h_length-lengthflippedx,1)];
            end
            y_est_1(ii,kk) = y_est_1(ii,kk) + sum(flippedx.*q_ES(:,jj));
            
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
plot(y_OS)
title('y_{sim1}')
%xlim([200 600])
% figure(14)
% plot(y_sim_2)
% title('y_{sim2}')
% %xlim([200 600])