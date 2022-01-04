hold on;
PTT_sum=0;
% cat doan tin hieu cua ECG sao cho bang voi PPG
new_ecg_signal=zeros(length(ppg_signal),1);
for b=1:length(ppg_signal)
    new_ecg_signal(b)=ecg_signal(b);
end
R_peak(R_peak>b)=[];
plot(new_ecg_signal);
plot(R_peak,new_ecg_signal(R_peak),'ro','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',6);% R_peak la vi tri con ecg_signal(R_peak) la gia tri
% nang cao muc tin hieu PPG
new_ppg_signal=zeros(length(ppg_signal),1);
for a=1:length(ppg_signal)
    new_ppg_signal(a)=ppg_signal(a)+600;
end
% ve tin hieu PPG
plot(new_ppg_signal);
plot(systole_peak,new_ppg_signal(systole_peak),'ro','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',6);
xlabel('Sample');
title('Systole points-PPG signal');
PTT_array=zeros(length(systole_peak),1);
for a=1+length(systole_peak)-length(R_peak):length(systole_peak)
    PTT_array(a)=(systole_peak(a)-R_peak(a-1))/fs;
    PTT_sum=PTT_sum+PTT_array(a);
end
PTT_average=PTT_sum/a;


