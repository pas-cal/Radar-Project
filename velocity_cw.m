function [] = velocity_cw(L, Tp, MS, norm)

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

close
clear 
clc
tic
%% Initialisations

[y, fs] = audioread('Ani Run Offset Start 1.wav');
c = 2.997e8;        % Speed of light
fc = 2.43e9;        % Carrier Central frequency
Tp = 0.1;           % DFT Pulse length Tp
N = fs*Tp;          % Number of samples in each pulse

v = -y(:,1);         % Negative Amplifier Gain
N_t = length(v);    % Number of Samples

%% Redundant initialisations

L = 4;
MS = 1;
norm = 2;

%% Raw Data Plot
figure
hold on;
plot(-y(:,1),'b--')
%plot(y(:,2),'r')
xlim([0.25*10^5 3.5*10^5])
ylabel('Amplitude')
xlabel('Data Samples')

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
%% Plot

vmax = c*fs/(4*fc);           % Velocity Resolution, fd_max = fs/2

% Axes resolution
time_array = linspace(0,Tp*count,count);
v_var = linspace(0,vmax,length(v_f_half));

figure
switch(norm)
    case 1
        imagesc(v_var,time_array,v_f_half,[-45,0])
    case 2
        imagesc(v_var,time_array,v_f_half,[-15,0])
end
xlabel('Velocity(m/s)')
xlim([0 15])
ylabel('Time(s)')
colorbar
fc = fc/10^9;
if MS == 0
    title('Continuous Wave Radar Velocity without MS Clutter Rejection')
else 
    title('Continuous Wave Radar Velocity with MS Clutter Rejection')
end
if norm == 1
    subtitle(['T_{p} =',num2str(Tp),'s; f_{c}',num2str(fc), 'GHz; ',num2str(L), '*N padding with Normalisation Method 1'])
else
    if norm == 2
        subtitle(['T_{p} =',num2str(Tp),'s; f_{c}',num2str(fc), 'GHz; ',num2str(L), '*N padding with Normalisation Method 2'])
    else
        subtitle(['T_{p} =',num2str(Tp),'s; f_{c}',num2str(fc), 'GHz; ',num2str(L), '*N padding without Normalisation'])
    end
end

v_f_half = v_f_half';
[~,M] = max(v_f_half);
v_space = v_var(M);
time_array = time_array(1:length(v_space));

v_scat_1 = zeros(1,size(v_f_half,2));
v_scat_2 = zeros(1,size(v_f_half,2));

for i = 1:size(v_f_half,2)
    [~,loc,w,p] = findpeaks(v_f_half(:,i));
    for j = 1:size(loc)-1
        maxi = j; 
        for k = j+1:size(loc)
            if v_f_half(loc(k),i) > v_f_half(loc(maxi),i)
                maxi = k;
            end
        end
        if maxi ~= j
            temp = loc(maxi);
            loc(maxi) = loc(j);
            loc(j)=temp;
            temp = w(maxi);
            w(maxi) = w(j);
            w(j)=temp;
            temp = p(maxi);
            p(maxi) = p(j);
            p(j)=temp;
        end
    end
    maxi = 2;
    for j =3:5 
        if w(j)*p(j) > w(maxi)*p(maxi)
            maxi = j;
        end
    end
    v_scat_1(1,i) = v_var(loc(1));
    v_scat_2(1,i) = v_var(loc(maxi));
end

windowSize = 10; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;

v_space = filter(b,a,v_space);
v_scat_1 = filter(b,a,v_scat_1);
v_scat_2 = filter(b,a,v_scat_2);


figure(3)
plot(time_array,v_space)
xlabel("Time [s]")
ylabel("Velocity (m/s)")
ylim([0 10])
title("Velocity vs Time Graph")


figure(4)
plot(time_array,v_scat_1,time_array,v_scat_2);
xlabel("Time [s]")
ylabel("Velocity (m/s)")
ylim([0 10])
title("Velocity vs Time Graph")
legend('Stongest', 'Second Strongest')
hold off
toc