% ecg =VarName1;
ecg=load('ecg _mr Tai.txt');
fs=1000;
Lfft=1024;
N1=100;
h1_ecg=fir1(N1,95*2/fs,'low'); % fir 100Hz
N2=500;
h2_ecg=fir1(N2,[52*2/fs 62*2/fs],'stop');% fir 50Hz
N3=100;
h3_ecg=fir1(N3,0.5*2/fs,'high');

L_ecg=length(ecg);
y1_ecg=conv(h1_ecg,ecg);
L1_ecg=length(y1_ecg);
y1_ecg=y1_ecg((N1+1)/2:L1_ecg-(N1+1)/2); %remove delay
y2_ecg=conv(h2_ecg,y1_ecg);
L2_ecg=length(y2_ecg);
y2_ecg=y2_ecg((N2+1)/2:L2_ecg-(N2+1)/2); %remove delay
y3_ecg=conv(y2_ecg,h3_ecg,'same');

figure;
subplot(211);plot(y2_ecg);
xlabel('Sample');
title('ECG before filter');
% subplot(212);plot(y2_ecg); 
subplot(212);plot(y3_ecg); 
xlabel('Sample');
title('ECG after filter');

figure;
subplot(211);plot(linspace(0,fs,Lfft),abs(fft(ecg,Lfft)));
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Before filter');
subplot(212);plot(linspace(0,fs,Lfft),abs(fft(y3_ecg,Lfft)));
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('After filter');

%===================R_detection on ECG===================%
ecg_signal=y3_ecg;
slp=slope(ecg_signal,12);% xet do doc cua tin hieu tham so dau vao: tin hieu va buoc nhay cua doc

high=high_detection(slp,200); % xac dinh diem dinh cua moi doc

min_high=10000;
for i=1:length(high)-3;
    high_point=slp(high(i));
    for j=i:i+3
        if slp(high(j))>high_point
            high_point=slp(high(j));
        end
    end
    if(high_point<min_high)
        min_high=high_point;
    end

end

threshold1=2*min_high;
noise_peak=high(slp(high)>threshold1);
high(slp(high)>threshold1)=[];

threshold2=0.76*(max(slp(high))+min(slp(high)))/2;
high(slp(high)<threshold2)=[];

peak=zeros(length(high),1);
n=75;
for i=1:length(high)
    max_point=ecg_signal(high(i));
    max_index=0;
    for j=-n:n
        if ecg_signal(high(i)+j)> max_point
            max_point=ecg_signal(high(i)+j);
            max_index=high(i)+j;
        end
    end
    peak(i)= max_index;
end
peak(peak==0)=[];

slp1=zeros(length(ecg_signal),1);
n1=30;
for i=1+n1:length(ecg_signal)
    if ecg_signal(i)> ecg_signal(i-n1)
        slp1(i)=ecg_signal(i)-ecg_signal(i-n1);
    else
        slp1(i)=0;
    end
end

slp1=slp1.^2;

threshold3=0.4*(min(slp1(peak))+ max(slp1(peak)))/2;
peak(slp1(peak)<threshold3)=[];
slp1(slp1<threshold3)=0;

R_peak=peak;

figure;
plot(ecg_signal);hold on;
plot(R_peak,ecg_signal(R_peak),'ro','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',6);
xlabel('Sample');
title('R points-ECG signal');

%===================Heart rate estimation===================%

HR=zeros(length(R_peak),1);
for i=2:length(R_peak)              % HR(1)=0
    j=find(noise_peak < R_peak(i),1,'last');
        if noise_peak(j)>R_peak(i-1)
            HR(i)=0;
        else
            HR(i)=60*fs/(R_peak(i)-R_peak(i-1));
        end
end

figure;
stem(HR,'ro','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',6);
title('Heart Rate');
