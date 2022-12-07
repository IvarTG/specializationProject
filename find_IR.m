function [ir] = find_IR(r_ES,n_meas,n_source,samplemethod,cair,fs)
%Finds the IR signal given the inputs and outputs of the system
% Discretization error for impulse responses

% samplemethod = 1 -> round at low fs
% samplemethod = 2 -> round at high fs and decimate
% samplemethod = 3 -> round at low fs but split pulse between two samples

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

end

