clear all;
close all;
clc;
avg_FC=zeros(50,68,68);
avg_ranFC=zeros(50,68,68);
sum_FC=zeros(68,68);
ransum_FC=zeros(68,68);
e=zeros(1,50);
e_ran=zeros(1,50);
matfiles = dir('/home/shruti/Desktop/Himadri/Shruti/IITG-final-work/Old-Controls/*_SC.mat');
matfiles=[matfiles; dir('/home/shruti/Desktop/Himadri/Shruti/IITG-final-work/Old-Controls/*_SC_new.mat')];
k=0;
ran_net=zeros(1,68,68);
for i=0.01:0.01:0.5
    k=k+1;
   for j=1:length(matfiles)  
    name=strcat('/home/shruti/Desktop/Himadri/Shruti/IITG-final-work/Old-Controls/',matfiles(j).name);
    load(name,'SC_cap_agg_bwflav1_norm');
    disp(matfiles(j).name);
    thresh_FC=threshold_proportional(SC_cap_agg_bwflav1_norm,i);
    adjFC=weight_conversion(thresh_FC,'binarize');
    sum_FC=sum_FC+adjFC;
    
    for l=1:1
    [R,~]=randmio_und(adjFC, 10);
    ran_net(l,:,:)=R;
    end
    r=randi(1,1,1);
    ransum_FC=ransum_FC+reshape(ran_net(r,:,:),68,68);
    clearvars -except sum_FC avg_FC matfiles i k j e thresh ran_net ransum_FC avg_ranFC e_ran
   end
    avg_FC(k,:,:)=sum_FC/length(matfiles);
    e(k)=entropy(avg_FC(k,:,:));
    
    avg_ranFC(k,:,:)=ransum_FC/length(matfiles);
    e_ran(k)=entropy(avg_ranFC(k,:,:));
end
   plot(1:50,e);
   hold on;
   plot(1:50,e_ran,'o--');
 thresh_diff=zeros(1,50); 
 for k=1:50 
     thresh_diff(k) = e(k)-e_ran(k); 
 end
 [~,thresh]=min(thresh_diff);
 thresh=thresh/100; %% 0.075                 
