function [time_space, v_space] = velocity_cw(L, Tp, MS, norm, y, fs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @ NAME: Process Velocity Data for CW radar
% @ INPUT:  L    ----- Pad Depth
%           Tp   ----- Pulse Time
%           MS   ----- MS Clutter Rejection Option [0 or 1]
%           norm  ----- Normalisation option - 1 or 2
% @ OUTPUT: -
% @ COMMENT: Normalisation 1 is subtracting max of matrix from the whole
% matrix, normalisation 2 is max of each row from that row.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialisations

c = 2.997e8;        % Speed of light
fc = 2.43e9;        % Carrier Central frequency
N = fs*Tp;          % Number of samples in each pulse

v = -y(:,1);         % Negative Amplifier Gain
N_t = length(v);    % Number of Samples

%% Raw Data Plot
% figure
% hold on;
% plot(-y(:,1),'b--')
% plot(y(:,2),'r')
% xlim([0.25*10^5 3.5*10^5])
% ylabel('Amplitude')
% xlabel('Data Samples')

%% Initialising the FFT Matrix in Time-domain

count = round(N_t/N);
v_t = y(1:(count-1)*(N),1);
v_t = reshape(v_t,[N,count-1])';


%% MS Clutter Rejection

switch(MS)
    case 1
         v_t = v_t-mean(mean(v_t)); % MS clutter removal
end

%% Zero Padding and frequency domain matrix

v_f = 20*log10(abs(fft(v_t,L*N,2)));

v_f_half = v_f(:,1:length(v_f)/2);      % Only lower half of the matrix is necesssary (frequency domain)

%% Normalisation

switch(norm)
    case 1
        v_f_half = v_f_half - max(max(v_f_half));
    case 2
        v_f_half = v_f_half - max(v_f_half,[],2);
end

time_space = linspace(0,Tp*count,count);

vmax = c*fs/(4*fc);  
v_var = linspace(0,vmax,length(v_f_half));


[~,M] = max(v_f_half');
v_space = v_var(M);
time_space = time_space(1:length(v_space));



