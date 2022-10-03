function [ir,fs] = readetxfile(filename,numberofsampelstoread);
% This function reads an EASERA etx file and returns the impulse
% response(s) stored therein and the sampling frequency. 
%
% Input parameters
%   filename            The name of the etx file, including a path so that
%                       Matlab can find the file.
%   numberofsampelstoread   (optional) If this value is specified, only
%                       this many sampels will be read for the impulse
%                       response(s) in the file.
%
% Output parameters
%   ir                  A matrix, size [nsampels,nchannels], with the
%                       impulse response(s) in the file.
%   fs                  The sampling frequency, in Hz.
%
% peter.svensson@ntnu.no 15 feb. 2022
%
% [ir,fs] = readetxfile(filename,numberofsampelstoread);

% 15 feb. 2022 First version

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The .etx files have 22 lines at the start which should be skipped
linestoskip = 22;

% .... but the .etx files have the the number of channels on line 7 
linewithnchannels = 7;

% .... and the sampling frequency on line 8
linewithfs = 8;

% .... and the number of samples on line 20
linewithnsamples = 20;

% We read just those lines first to get this data

fid = fopen(filename,'r');
for jj = 1:linewithnchannels
    textlinenchannels = fgetl(fid);
end
for jj = 1:linewithfs-linewithnchannels
    textlinefs = fgetl(fid);
end
for jj = 1:linewithnsamples-linewithfs
    textlinensamples = fgetl(fid);
end
fclose(fid);

iv = strfind(textlinenchannels,'Channels');
numstring = textlinenchannels(iv+9:end);
nchannels = str2num(numstring);

iv = strfind(textlinefs,'[Hz]');
numstring = textlinefs(iv+4:end);
fs = str2num(numstring);

iv = strfind(textlinensamples,'TimeSamples');
numstring = textlinensamples(iv+11:end);
nsamplesinfile = str2num(numstring);

if nargin == 1
    numberofsampelstoread = nsamplesinfile;
end
if numberofsampelstoread > nsamplesinfile
    error('ERROR: The number of samples in the file is smaller than what you asked for')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

endlinetoread = linestoskip + 1 + numberofsampelstoread - 1;

% The first impulse response will be stored in column B = column 2.
startcolumn  ='B';
endcolumn = setstr(startcolumn + nchannels - 1);

% linereadtostring specifies which rows and columns should be read in the
% text file.
linereadstring = [startcolumn,int2str(linestoskip+1),':',endcolumn,int2str(endlinetoread)];

ir = readmatrix(filename,'FileType','text','Range',linereadstring);
