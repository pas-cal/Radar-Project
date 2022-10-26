
clc, clear;


%% Initialisation
y = [];
time_space = 0;
range_space = 0;

Fs = 44100 ; 
bit_depth = 16 ; 
channels = 2 ; 
ID = -1;             % Have to check which channel it is
input_stream = audiorecorder(Fs,bit_depth,channels,ID);

windowSize = 5; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;

T_s = 0.2;
x = 0;

flag = 1;


%% Realtime Loop

while(flag == 1)

    x = x + T_s; 
    recordblocking(input_stream,T_s);
    audioData = input_stream.getaudiodata;
    y = [y;audioData];

    [range_space, time_space] = Range_Data(4, 20*10^-3, 1, 3, y, Fs);

    figure(1)
    plot(time_space,range_space,'o')
    xlabel("Time [s]")
    ylabel("Range (m)")
    ylim([0 40])
   
    
    ButtonHandle = uicontrol('Style', 'PushButton', 'String', 'Stop Measurement','Callback', 'delete(gcbf)');
    uicontrol(ButtonHandle)
      if ~ishghandle(ButtonHandle)
          flag = 0;
          disp('Loop stopped by user');
          break;
      end
end

figure(2)
plot(y(:,1))
title('Sampled Signal')
xlabel('Sample Number')

%% Plot All Sampled Data
figure(3)
plot(time_space,range_space)
xlabel("Time [s]")
ylabel("Range (m)")
ylim([0 40])
title("Final Recorded Signal")

% figure(4)
% imagesc(range_space,time_space,chirp_freq_half,[-50,0])
% xlabel('Range (m)')
% ylabel('Time (s)')
% xlim([0 100]);