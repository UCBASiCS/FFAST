%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab script to demonstrate the workings of FFAST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;

%% CONFIGURATIONS
% coprime factors
coPrimeFactors = [49 50 51];
% Fourier sparsity
k = 5;
% noisy or not
isNoisy = 1;                    %% <= CHANGE HERE IF NOISY SIGNAL IS WANTED
% if noisy is chosen, the SNR
SNRDB = 10;
% infile
infile = 'infile.txt';
% true sinusiods
tsfile = 'truesinusoids.txt';

% resulting signal length from coprime factors
n = prod(coPrimeFactors);

% the fundamental frequency
omega = 2*pi/n;

% the frequency indices
fqsIdx = randsample(n,k);
fqsIdx = sort(fqsIdx);
% the frequency magnitudes
amps = ones(k,1) + rand(k,1);

fqs = fqsIdx*omega;

% generate signal vector
disp('generating signal...');
signal = amps'*exp(1i*fqs*(0:n-1));

x2Min = min(amps.^2);

% if isNoisy add noise to the signal
if isNoisy
    disp('adding noise...');
    SNR = 10^(SNRDB/10);
    sigma2 = x2Min/(SNR*n);
    noise = sqrt( sigma2/2 ) *( randn(1,n) + 1i*randn(1,n) );
    signal = signal + noise;
end

% write the signal file
disp('writing time domain signal to file...');
fid = fopen(infile, 'w');
for i = 1:n
    fprintf(fid, '(%f,%f)\n', real(signal(i)), imag(signal(i)));
end
fclose(fid);

% write the true values
disp('writing the frequency domain signal to file...');
fid = fopen(tsfile, 'w');
fprintf(fid,'(frequency bin, amplitude)\n');
fprintf('frequency bin, magnitude and phase\n');
for i = 1:k
    fprintf(fid, '(%d,%f)\n', fqsIdx(i), amps(i));
    fprintf('%8d %8.3f %8.3f\n', fqsIdx(i), abs(amps(i)), angle(amps(i)) );
end
fclose(fid);

%% RUN FFAST
disp('running FFAST...');
if isNoisy
    system(sprintf('../exe/ffast -f infile.txt -k %d -s %d', k, SNRDB))
else
    system(sprintf('../exe/ffast -f infile.txt -k %d',k))
end
%% DISPLAY TRUE SPECTRUM AND FFAST OUTPUT
disp('parsing FFAST output file...');
recSig = parse_output();

stem(fqsIdx, amps, 'b');
hold on
if size(recSig,1) >=1
    plot(recSig(:,1), abs(recSig(:,2))/sqrt(n),'rx');
else
    plot(0, 0,'rx');
end
hold on
    
if isNoisy
    title(sprintf('FFAST recovery, n = %d, k = %d, SNRdB = %.2f ', n, k, SNRDB));
else
    title(sprintf('FFAST recovery, n = %d, k = %d', n, k));
end

legend('true spectrum','FFAST output');
legend('Location','southoutside','Orientation','horizontal');
legend('boxoff');
xlabel('frequency');
ylabel('magnitude');

disp('== FFAST output ==');
fprintf('frequency, magnitude and phase\n');
recSig = sortrows(recSig);
if size(recSig,1) >=1
    for i = 1:size(recSig,1)
        fprintf('%8d %8.3f %8.3f\n',recSig(i,1),abs(recSig(i,2)/sqrt(n)),angle(recSig(i,2)));
    end
end