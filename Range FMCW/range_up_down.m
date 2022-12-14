function [] = range_up_down(L, Tp, MS, MTI)

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

[y,fs] = audioread('1 person multi speed walk test 1.wav');
N_t = length(y);                % Total number of samples in positive chrip read file
data = -y(:,1);                 % Negative of Raw data because of the use of an Inverting Amplifier
sync = -y(:,2);                 % data array - backscatter information; sync array - sync from modulator

f_start = 2.405*10^9;           % Start frequency
f_stop = 2.495*10^9;            % Stop frequency
c = 2.997*10^8;                 % Speed of Light (m/s)

band = f_stop-f_start;          % Bandwidth of the signal
range_res = c/(2*band);         % Range Resolution

f_c = band/2;
lambda = c/f_c;

Tp = 20e-3;
N=Tp*fs;                        % Number of samples in each pulse segment
T_max = N_t/fs;                 

%% Variable declarations for non function use

MS = 1; 
MTI = 3;
L = 4;

%% Raw Data Plot
figure (1)
hold on;
xlim([0.25*10^5 3.5*10^5])
ylabel('Amplitude')
xlabel('Data Samples')
subplot(2,1,1)
plot(-y(:,1),'b--')
title('Data Captured')
subplot(2,1,2)
plot(y(:,2),'r')
title('Sync')
hold off

%% Up-chirp Processing

chirp_logic = (sync > 0.1);       % Logical Array of each sample when modulator sync amplitude is positive
x = 0;
for i = 2:(size(chirp_logic)-N)
    if chirp_logic(i) == 1 && chirp_logic(i-1) == 0 && chirp_logic(i-2) == 0      % To make sure we get a proper chirp in the sample and not a random segment
        x = x + 1;
        chirp_time(x,:) = data(i:i+N-1);
    end
end
cont = reshape(chirp_time,1,[]);

figure (2)
hold on;
xlim([0.25*10^5 3.5*10^5])
ylabel('Amplitude')
xlabel('Data Samples')
plot(cont,'b--')
title('Parsed Data Captured')
plot(y(:,2),'r')
legend('Parsed Data','Sync')
hold off

%% Up and Down Chirp

up_column = 0; up_row = 1; down_row = 1; down_column = 0;
for i = 1:(length(data)-1)
    if chirp_logic(i) == 1
        up_column = up_column + 1;
        updata_mat(up_row,up_column) = data(i);
        if chirp_logic(i+1) ~= 1
            up_row = up_row + 1;
            up_column = 0;
        end
    elseif chirp_logic(i) == 0
        down_column = down_column + 1;
        downdata_mat(down_row,down_column) = data(i);
        if chirp_logic(i+1) ~= 0
            down_row = down_row + 1;
            down_column = 0;
        end
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
        for i = 1:N
            Ms_mean = mean(updata_mat(:,i));
            updata_mat(:,i) = updata_mat(:,i) - Ms_mean;
        end
        for i = 1:N
            Ms_mean = mean(downdata_mat(:,i));
            downdata_mat(:,i) = downdata_mat(:,i) - Ms_mean;
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
        updata_mat = updata_mat(3:size(updata_mat,1),:)-updata_mat(2:size(updata_mat,1)-1,:)-(updata_mat(2:size(updata_mat,1)-1,:)-updata_mat(1:size(updata_mat,1)-2,:));
        downdata_mat = downdata_mat(3:size(downdata_mat,1),:)-downdata_mat(2:size(downdata_mat,1)-1,:)-(downdata_mat(2:size(downdata_mat,1)-1,:)-downdata_mat(1:size(downdata_mat,1)-2,:));
        x = x - 2;
end

%% Performing zero padding and IFFT on time domain matrix
for i = 1:x
    chirp_freq(i,:) = 20*log10(abs(ifft(chirp_time(i,:),L*N)));        % Calculate IFFT using zero padding
    updata_freq(i,:) = 20*log10(abs(ifft(updata_mat(i,:),L*N)));  
    downdata_freq(i,:) = 20*log10(abs(ifft(downdata_mat(i,:),L*N)));
end   

fmax = length(chirp_freq)/2;                                        % Nyquist
chirp_freq_half = chirp_freq(:,1:fmax);                             % Generate the first half part
chirp_freq_half = chirp_freq_half - max(max(chirp_freq_half));      % Normalization

fmax = length(updata_freq)/2;                                        % Nyquist
updata_freq_half = updata_freq(:,1:fmax);                             % Generate the first half part
updata_freq_half = updata_freq_half - max(max(updata_freq_half));      % Normalization

fmax = length(downdata_freq)/2;                                        % Nyquist
downdata_freq_half = downdata_freq(:,1:fmax);                             % Generate the first half part
downdata_freq_half = downdata_freq_half - max(max(downdata_freq_half));      % Normalization

range_t = linspace(0,(range_res*N/2),L*N);            % Create frequency array
time_array = linspace(0,T_max,x);                   % Time Array

[n,k]= size(downdata_freq_half);
for i=2:n
    [~,n1] = max(downdata_freq_half(i,:));
    [~,n2] = max(updata_freq_half(i,:));
    if(abs(n1-n2)>5)
        index(i)=index(i-1);
    else 
        index(i)=n1-n2;
    end
    vel_cw(i)= (n1-n2)*fs/k;
end

velocity_array = lambda*vel_cw/400;
%% Plot and Misc
figure(3)
imagesc(range_t/2,time_array,chirp_freq_half,[-50,0])
xlabel('Range (m)')
ylabel('Time (s)')
xlim([0 50]);
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

%% Strongest scatterers

chirp_freq_half = chirp_freq_half';
[~,M] = max(chirp_freq_half);
range = range_t(M);
range = range(1:length(time_array));

r_scat_1 = zeros(1,size(chirp_freq_half,2));
r_scat_2 = zeros(1,size(chirp_freq_half,2));

for i = 1:size(chirp_freq_half,2)
    [~,loc,w,p] = findpeaks(chirp_freq_half(:,i));
    for j = 1:size(loc)-1
        maxi = j; 
        for k_u = j+1:size(loc)
            if chirp_freq_half(loc(k_u),i) > chirp_freq_half(loc(maxi),i)
                maxi = k_u;
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
    for j =3:13 
        if w(j)*p(j) > w(maxi)*p(maxi)
            maxi = j;
        end
    end
    r_scat_1(1,i) = range_t(loc(1));
    r_scat_2(1,i) = range_t(loc(maxi));
end

windowSize = 10; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;

range = filter(b,a,range);
r_scat_1 = filter(b,a,r_scat_1);
r_scat_2 = filter(b,a,r_scat_2);



windowSize = 50; 
b = (1/windowSize)*ones(1,windowSize);
velocity_array = filter(b,a,velocity_array);
for i=2:length(velocity_array)
    if (abs(velocity_array(i)-velocity_array(i-1))>6) % Condition max velocity can not be more than 6 m/s for normal human
        velocity_array(i)=velocity_array(i-1);
    end
end


hold off
figure(4)
plot(time_array,range)
xlabel("Time [s]")
ylabel("Range (m)")
ylim([0 40])
xlim([0 25])
title("Range vs Time Graph")

figure(5)
plot(time_array,r_scat_1,time_array,r_scat_2);
xlabel("Time [s]")
ylabel("Range (m)")
ylim([0 40])
title("Range vs Time Graph")
legend('Stongest', 'Second Strongest')
hold off
% 
% figure(6)
% plot(time_array_v,velocity);
% ylabel('Velocity(m/s)');
% xlim([0 25]);
% xlabel('Time(s)')
% title("Velocity (finite difference) vs Time Graph")

%% Up and down chirp plot
figure(7)
subplot(1,2,1);
imagesc(range_t/2,time_array,updata_freq_half,[-50,0])
xlabel('Range(m)')
ylabel('Time(s)')
title('Upchirp Range Data');
subtitle(['T_{p} =',num2str(Tp),'ms; f_{start}',num2str(f_start), 'GHz; f_{stop}',num2str(f_stop),'GHz; ',num2str(L), '*N padding'])
xlim([0 40])
colorbar

subplot(1,2,2);
imagesc(range_t/2,time_array,downdata_freq_half,[-50,0])
xlabel('Range(m)')
ylabel('Time(s)')
title('Downchirp Range Data');
xlim([0 40])
colorbar
subtitle(['T_{p} =',num2str(Tp),'ms; f_{start}',num2str(f_start), 'GHz; f_{stop}',num2str(f_stop),'GHz; ',num2str(L), '*N padding'])

figure(8)
plot(time_array,velocity_array);
ylabel('Velocity(m/s)');
xlim([0 25]);
xlabel('Time(s)')
hold off
title('Velocity From Simultaneous Processing')
toc