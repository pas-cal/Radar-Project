function [] = range_cw(L, Tp, MS, MTI)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @ NAME: Process Range Data for FMCW radar
% @ INPUT:  L    ----- Pad Depth
%           Tp   ----- Pulse Time
%           MS   ----- MS Clutter Rejection Option [0 or 1]
%           MTI  ----- MTI Clutter Rejection Option [0 - off or 2 or 3]
% @ OUTPUT: -
% @ COMMENT: -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear
clc

tic
%% Initialisation Phase

[y,fs] = audioread('Range_Test_File.mp4');
N_t = length(y);                % Total number of samples in positive chrip read file
data = -y(:,1);                 % Negative of Raw data because of the use of an Inverting Amplifier
sync = -y(:,2);                 % data array - backscatter information; sync array - sync from modulator

f_start = 2.35*10^9;           % Start frequency
f_stop = 2.55*10^9;            % Stop frequency
c = 2.997*10^8;                 % Speed of Light (m/s)

band = f_stop-f_start;          % Bandwidth of the signal
range_res = c/(2*band);         % Range Resolution

Tp = 20e-3;
N=Tp*fs;                        % Number of samples in each pulse segment
T_max = N_t/fs;                 

%% Variable declarations for non function use

MS = 0;
MTI = 2;
L = 4;
%% Up-chirp Processing

chirp_logic = (sync > 0);       % Logical Array of each sample when modulator sync amplitude is positive
x = 0;
for i = 1:(size(chirp_logic)-N)
    if chirp_logic(i) == 1 && chirp_logic(i-1) == 0 && chirp_logic(i-2) == 0      % To make sure we get a proper chirp in the sample and not a random segment
        x = x + 1;
        chirp_time(x,:) = data(i:i+N-1);
        time(x) = i/fs;
    end
end
%% MS Clutter Reduction Step
if MTI == 0
    switch(MS)
        case 1
        for i = 1:N
            Ms_mean = mean(chirp_time(:,i));                  % Mean accross each column 
            chirp_time(:,i) = chirp_time(:,i) - Ms_mean;        % MS clutter Reduction
        end
    end
end

%% MTI Step
switch(MTI)
    case 2
        chirp_time = chirp_time(2:size(chirp_time,1),:)-chirp_time(1:size(chirp_time,1)-1,:);
        x = x - 1;
    case 3
        chirp_time = chirp_time(3:size(chirp_time,1),:)-chirp_time(2:size(chirp_time,1)-1,:)-(chirp_time(2:size(chirp_time,1)-1,:)-chirp_time(1:size(chirp_time,1)-2,:));
        x = x - 2;
end

%% Performing zero padding and IFFT on time domain matrix
for i = 1:x
    chirp_freq(i,:) = 20*log10(abs(ifft(chirp_time(i,:),L*N)));        % Calculate IFFT using zero padding
end   
fmax = length(chirp_freq)/2;                                         % Nyquist
chirp_freq_half = chirp_freq(:,1:fmax);                           % Generate the first half part
chirp_freq_half = chirp_freq_half - max(max(chirp_freq_half)); % Normalization
range = linspace(0,(range_res*N/2),L*N);          % Create frequency array
time_array = linspace(0,T_max,x);    % Time Array

figure
imagesc(range,time_array,chirp_freq_half,[-50,0])
xlabel('Range (m)')
ylabel('Time (s)')
xlim([0 100]);
colorbar
if MS == 0 && MTI == 0
    title('Continuous Wave Radar Range Measurement without Clutter Rejection')
else 
    if MTI == 0
        title('Continuous Wave Radar Range Measurement with MS Clutter Rejection')
    end
end
if MTI == 2
    title('Continuous Wave Radar Range Measurement with 2-pulse MTI Clutter Rejection')
else
    if MTI == 3
        title('Continuous Wave Radar Range Measurement with 3-pulse MTI Clutter Rejection')
    end
end
Tp = Tp*10^3;
f_start = f_start/10^9;
f_stop = f_stop/10^9;
subtitle(['T_{p} =',num2str(Tp),'ms; f_{start}',num2str(f_start), 'GHz; f_{stop}',num2str(f_stop),'GHz; ',num2str(L), '*N padding'])

hold on
toc
