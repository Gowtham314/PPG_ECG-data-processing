close all;
clear all;
clc

fs=1000;
index=31;
 % load du lieu
ecg=load('D:\Matlab\ppg-mr the anh.txt');
ppg=load('D:\Matlab\ecg-mr the anh.txt');

%===================ECG filter===================%
Lfft=4096;
N1=100;
h1_ecg=fir1(N1,95*2/fs,'low'); % fir 100Hz
N2=500;
h2_ecg=fir1(N2,[38*2/fs 45*2/fs],'stop');% fir 50Hz
N3=100;
h3_ecg=fir1(N3,0.5*2/fs,'high');

 % Loc va remove delay
y1_ecg=conv(ecg(1:length(ecg)),h1_ecg,'same');
y2_ecg=conv(y1_ecg(1:length(ecg)),h2_ecg,'same');
y3_ecg=conv(y2_ecg(1:length(ecg)-50),h3_ecg,'same');

%======== ve do thi tin hieu dien tim truoc va sau khi loc ======%
figure;
subplot(211);plot(ecg);
xlabel('Sample');          
title('ECG before filter');
subplot(212);plot(y3_ecg); 
xlabel('Sample');
title('ECG after filter');
%======= ve pho bien do tan so cua tin hieu ECG truoc va sau khi loc=======%
% figure;
% subplot(211);plot(linspace(0,fs,Lfft),abs(fft(ecg,Lfft)));      % 
% xlabel('Frequency (Hz)');
% ylabel('Amplitude');
% title('Before filter');
% subplot(212);plot(linspace(0,fs,Lfft),abs(fft(y3_ecg,Lfft)));
% xlabel('Frequency (Hz)');
% ylabel('Amplitude');
% title('After filter');

%===================R_detection on ECG===================%
ecg_signal=y3_ecg;  % 
slp=slope(ecg_signal,12);% tinh do doc cua ECG voi buoc nhay la 12 de lam giam truong hop cua nhieu
high=high_detection(slp,200); % Tim vi tri cua dinh doc voi buoc nhay la 200 (chon 150ms de nho hon 250ms (R-R)) 

min_high=1000;
%============
%======= Tao 1 cua so 3 gia tri de tim gia tri lon nhat cua cua so vao bien high_point va cho dich nhu thuat toan high_detection===========
for i=1:length(high)-3;
    high_point=slp(high(i)); 
    for j=i:i+3              
        if slp(high(j))>high_point   %xet xem trong trong 3 phan tu tiep theo phan tu nao lon hon gia tri dinh
          high_point=slp(high(j));  % thi gan diem do cho high_point
        end
    end
    %==== cap nhap gia tri thread_hold moi la high_port khi highport<=============
    if(high_point<min_high)         
         min_high=high_point;      
    end
end

threshold1=2*min_high; % dat muc nguong toi thieu cho do doc
noise_peak=high(slp(high)>threshold1); % detect dinh doc nhieu cao bat thuong (do nhieu co)
high(slp(high)>threshold1)=[];  %  cac vi tri ma tai do dinh doc lon hon nguong nao do
% Loai bo cac thanh dinh doc duoi muc threadhold 2
threshold2=0.76*(max(slp(high))+min(slp(high)))/2; % day la gia tri trung binh
high(slp(high)<threshold2)=[];  % nhung dinh doc nao nho hon muc nguong thu 2 thi xoa no di
% luc nay high chua cac index cua cac dinh doc dat tieu chuan
peak=zeros(length(high),1);  % khoi tao index cua cac diem dinh co do dai bang do dai cua mang high(so dinh doc)
n=75;   % 
%============== tim vi tri dinh================

for i=1:length(high)
    max_point=ecg_signal(high(i));  
    max_index=0;  
    % chon n=75 vi dinh o lan can dinh doc ve phia truoc hoac sau 75 diem (can co bai bao cho so lieu nay) 
    % khoang cach giua hai dinh toi thieu phai dat 150 mau
    for j=-n:n
        if ecg_signal(high(i)+j)> max_point
            max_point=ecg_signal(high(i)+j);
            max_index=high(i)+j;  % tra ve vi tri dinh trong truong hop trong cua so ay co 1 dinh
        end
    end
    peak(i)= max_index;  % tra ve gia tri 0 neu vi tri dinh doc khong the la diem dinh
end
peak(peak==0)=[]; % loai bo cac phan tu 0 

slp1=zeros(length(ecg_signal),1);
n1=30;
% tim do doc cua tin hieu ecg voi buoc nhay 30
for i=1+n1:length(ecg_signal)
    if ecg_signal(i)> ecg_signal(i-n1)
        slp1(i)=ecg_signal(i)-ecg_signal(i-n1);
    else
        slp1(i)=0;
    end
end
% binh phuong gia tri cac phan tu cua mang do doc moi
slp1=slp1.^2;
% threashold 3 la gia tri do doc toi thieu ma tai do 1 dinh R phai dat duoc
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
HR_sum=0;
for i=2:length(R_peak)              % HR(1)=0
    j=find(noise_peak < R_peak(i),1,'last');
        if noise_peak(j)>R_peak(i-1)
            HR(i)=0;
        else
            HR(i)=60*fs/(R_peak(i)-R_peak(i-1));
        end
     HR_sum=HR_sum+HR(i);   
end
HR_average=HR_sum/(length(R_peak)-1);
% figure;
% stem(HR,'ro','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',6);
% title('Heart Rate');

%===================PPG filter===================%
N1=40;
h1_ppg=fir1(N1,20*2/fs,'low'); % fir 20Hz
N2=80;
h2_ppg=fir1(N2,0.2*2/fs,'high'); % fir 0.2Hz
y1_ppg=conv(ppg,h1_ppg,'same');
y2_ppg=conv(y1_ppg(1:length(ppg)-10),h2_ppg,'same');
% L_ppg=length(ppg);
% y1_ppg=conv(h1_ppg,ppg);
% L1_ppg=length(y1_ppg);
% y1_ppg=y1_ppg((N1+1)/2:L1_ppg-(N1+1)/2); %remove delay

% y2_ppg=conv(h2_ppg,y1_ppg);
% L2_ppg=length(y2_ppg);
% y2_ppg=y2_ppg((N2+1)/2:L2_ppg-(N2+1)/2); %remove delay



figure;
subplot(211);plot(ppg);
xlabel('Sample');
title('PPG before filter');
subplot(212);plot(y2_ppg);
xlabel('Sample');
title('PPG after filter');
% ylim([1.5 3.5]);
%=========Ve fft cua tin hieu ecg================
% figure;
% subplot(211);plot(linspace(0,fs,Lfft),abs(fft(ppg,Lfft)));
% xlabel('Frequency (Hz)');
% ylabel('Amplitude');
% title('Before filter');
% subplot(212);plot(linspace(0,fs,Lfft),abs(fft(y2_ppg,Lfft)));
% xlabel('Frequency (Hz)');
% ylabel('Amplitude');
% title('After filter');
%===================Max_slope_detection on PPG===================%
ppg_signal=y2_ppg;
slp=slope(ppg_signal,12); % tim doc cua ppg voi buoc nhay 12
% figure;
% plot(slp);

high=high_detection(slp,200);  % xac dinh dinh cua doc nhu tin hieu ecg
% figure;
% subplot(211);
% plot(slp);hold on;
% plot(high,slp(high),'ro','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',6);
% 
% subplot(212);
% plot(ppg_signal);hold on;
% plot(high,ppg_signal(high),'ro','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',6);
%=========== Xac dinh Threadhold1 de loai bo cac dinh doc lon qua muc do nhieu ============%
min_high=1000;
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
 threshold1=4*min_high;
 %=======Loai nhieu==========
 high(slp(high)>threshold1)=[]; %remove detected points resulting from moving noise
 %=======Ve hinh=========
 figure;
 plot(slp);hold on;
 plot(high,slp(high),'ro','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',6);


threshold2=0.7*(max(slp(high))+min(slp(high)))/2;
high(slp(high)<threshold2)=[];
% figure;
% plot(slp);hold on;
% plot(high,slp(high),'ro','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',6);
%=========== Ve tin hieu=============
%============ PPG_diff la mang vi tri cua suon xung ======
PPG_diff=high;     % vi tri dinh cua tin hieu doc la vi tri suon cua tins hieu PPG
% PPG cung voi cac diem suon
figure;
plot(ppg_signal); hold on;
plot(PPG_diff,ppg_signal(PPG_diff),'ro','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',6);
xlabel('Sample');
title('Max slope-PPG signal');


%===================PTT_estimation===================%
PTT=zeros(length(R_peak),1);
PTT_sum=0;
for i=1:length(R_peak)
    j=find(PPG_diff > R_peak(i),1);
    if j>0
        PTT(i)=(PPG_diff(j)-R_peak(i))/fs;
    else
        PTT(i)=0;                                     %Consider the leftmost R_peak lack of corresponding PPG_diff
    end
    PTT_sum=PTT_sum+PTT(i);
end
PTT_average=PTT_sum/length(R_peak);
%=========== Ve tin hieu ECG va PPG cung cac diem suon (PPG) va cac diem dinh R (ECG)=============%
figure;
plot(ecg_signal-1000);hold on;
plot(R_peak,ecg_signal(R_peak)-1000,'ro','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',6);
plot(ppg_signal+500);
stem(PPG_diff,ppg_signal(PPG_diff)+500,'ro','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',6);
xlabel('Sample');
title('ECG-PPG');

%=======Ve PTT==============
figure;
stem(PTT,'ro','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',6);
xlabel('Sample');
ylabel('Second(s)');
title('PTT');

