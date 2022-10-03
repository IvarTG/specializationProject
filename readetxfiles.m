% avstandsmalinger_time.m
%
% Analyze impulse responses measured along a straight line in front of the small
% tube loudspeaker.
%
% Measurements done by Olav Ellingsen 12 Nov. 2020
% Analysis by Peter Svensson 12 Nov. 2020

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Values set by the user

% The file names
nfiles = 43;
filestart = 'M000';
fileend = '_S01_R01.etx';

% The microphone positions
xmic = [2.7:-0.1:-1.5];

% Choose how many samples of the IRs should be saved in the mat file and in
% the allirs matrix.
irlengthtosave = 44100;

% Choose low-pass filter for the lp-filtered versions of the impulse
% responses. We want to low-pass filter so that we get a truly
% omni-directional radiation. Suggestion: lpfreq = 2000 and nfir = 30

% The cut-off frequency for the low-pass filter, and the number of fir
% filter coefficients
lpfreq = 20000;
% nfir = 30;
% nfir = 10;
nfir = 6;

cair = 343;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the files
% The first time the script is run (or if you have cleared the workspace), 
% the .etx files will be loaded. They will be saved in a file allirs.mat,
% and after that first run, this file (allirs.mat) will automatically be loaded.
   
if exist('loadorigfiles') == 0
    loadorigfiles = 1;
end

if loadorigfiles == 1
    for ii = 1:nfiles
       if ii < 10 
           filename = [filestart,int2str(ii),fileend];
       else
           filename = [filestart(1:end-1),int2str(ii),fileend];           
       end
       fid=fopen(filename);
        cdata1=textscan(fid,'%f%f','delimiter',',', 'HeaderLines', 22 );
        fclose(fid);    
        ir1 = cdata1{2};
        
        if ii == 1
            allirs = zeros(irlengthtosave,nfiles);  
            tvec = cdata1{1};
            fs = 1/(tvec(2)-tvec(1));
        end
        allirs(:,ii) = ir1(1:irlengthtosave);
    %    
    end
    
    save allirs.mat allirs fs
    clear ir1 tvec
    loadorigfiles = 0;
else
   load allirs.mat   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create low-pass filtered versions of the impulse responses

Blp = fir1(nfir,lpfreq/(fs/2));
% Blp = fir1(nfir,lpfreq*[1/sqrt(2) sqrt(2)]/(fs/2));

[nirlength,nirs] = size(allirs);
allirslp = zeros(nirlength+nfir,nirs);

for ii = 1:nirs
   allirslp(:,ii) = conv(allirs(:,ii),Blp); 
end
allirslp = allirslp(nfir/2+1:end,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot one example of broadband IR and lp-filtered IR

tvec = 1/fs*[0:size(allirs,1)-1];

figure(1)
clf(1)
iv = 1:300;
h = plot(tvec(iv)*1e3,allirs(iv,1),'-o',tvec(iv)*1e3,10*allirslp(iv,1),'-o');
seth
set(h(1),'LineWidth',2);
set(h(2),'LineWidth',2);
g = xlabel('Time   [ms]');
sg14
grid
g = ylabel('Impulse response   [-]');
sg14
g = legend('Broadband IR',['LP-filtered IR (',num2str(round(lpfreq/100)/10),'kHz)']);
sg14
g = title('Example of measured impulse responses');
sg14

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot all the lp-filtered IR

figure(4)
h = plot(tvec*1e3,allirslp(1:irlengthtosave,:));
seth
g = xlabel('Time   [ms]');
sg14
grid
g = ylabel('Impulse response   [-]');
sg14
g = title('All measured impulse responses, lp-filtered');
sg14

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find all the peak values of the lp-filtered IRs

nmic = size(allirs,2);

allpeakvals = zeros(1,nmic);

for ii = 1:nmic
    [npeaklp,ninit] = findinit(allirslp(:,ii).^2,0.001);
    allpeakvals(ii) = allirslp(npeaklp,ii);
end

xplot = [0:0.001:max(xmic)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make a 2nd-order polynomial fit to
% 1/IRpeak^2. 

% [polyfitparams_amplin,S_amplin] = polyfit(xmic,1./allpeakvals,1);
% smoothcurve_amplin = polyfitparams_amplin(2) + polyfitparams_amplin(1)*xplot;
% xs_lin = -polyfitparams_amplin(2)/polyfitparams_amplin(1);
% alpha = 0.05; % Significance level
% [yfit,delta_lin] = polyconf(polyfitparams_amplin,xs_lin,S_amplin,'alpha',alpha);
% xs_lin_predint = delta_lin/polyfitparams_amplin(1);

[polyfitparams_ampsquare,S_ampsquare] = polyfit(xmic,1./allpeakvals.^2,2);
smoothcurve_ampsquare = polyfitparams_ampsquare(3) + polyfitparams_ampsquare(2)*xplot + polyfitparams_ampsquare(1)*xplot.^2;
xs_square =  -polyfitparams_ampsquare(2)/polyfitparams_ampsquare(1)/2;
dist2_square = polyfitparams_ampsquare(3)/polyfitparams_ampsquare(1) - xs_square^2;
dist_square = sqrt(dist2_square);
ampsquaremin = dist2_square*polyfitparams_ampsquare(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the results for the fit to 1/IRpeak^2

figure(3)
h = plot(xmic,1./allpeakvals.^2,'o',xplot,smoothcurve_ampsquare,'--',...
    xs_square,ampsquaremin,'k*',xs_square,ampsquaremin,'k*');
seth
set(h(1),'LineWidth',2);
set(h(2),'LineWidth',2);
set(h(3),'LineWidth',2,'MarkerSize',12);
grid
g = xlabel(['xmic   [m]']);
sg14
g = ylabel(['1/IR peak amplitude squared   [-]']);
sg14
g = legend('Data from lp-filtered IR (2k)','Curve fitting',...
    ['Est. xs=',num2str(round(1e4*xs_square)/10),'mm'],...
    ['Est. yz-dist=',num2str(round(dist_square*1e4)/10),'mm']);
sg14
set(g,'Location','best')
g = title('Inverse of IR peak amplitude, squared');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

iv = 100:700;

plotscaling = 0.3*lpfreq/1e4;

plotshiftmatrix = [1:nfiles];
plotshiftmatrix = plotshiftmatrix(ones(length(iv),1),:);
figure(4)
plot(tvec(iv)*1e3,allirslp(iv,:)+plotscaling*plotshiftmatrix)
grid





