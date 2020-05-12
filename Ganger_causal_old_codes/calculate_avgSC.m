clear all;
close all;
clc;
%%% YOUNG CONTROLS %%%
matfiles = dir('/home/shruti/Desktop/Himadri/Shruti/IITG-final-work/Young-Controls/*_SC.mat'); 
matfiles=[matfiles;dir('/home/shruti/Desktop/Himadri/Shruti/IITG-final-work/Young-Controls/*_SC_new.mat')];
avg_SC=zeros(68,68);
for i=1:length(matfiles)  
    name=strcat('/home/shruti/Desktop/Himadri/Shruti/IITG-final-work/Young-Controls/',matfiles(i).name);
    load(name,'SC_cap_agg_bwflav1_norm');
    disp(matfiles(i).name);
    thresh_SC=threshold_proportional(SC_cap_agg_bwflav1_norm,0.11); % 30 percent of strongest weights preserved
    avg_SC=avg_SC+thresh_SC;
    save(name,'thresh_SC','-append');
    clearvars -except avg_SC matfiles i
end
avg_SC_young=avg_SC/length(matfiles);
save('avg_SC_young.mat','avg_SC_young');

%%% OLD CONTROLS %%%%

 matfiles = dir('/home/shruti/Desktop/Himadri/Shruti/IITG-final-work/Old-Controls/*_SC.mat'); 
 matfiles=[matfiles;dir('/home/shruti/Desktop/Himadri/Shruti/IITG-final-work/Old-Controls/*_SC_new.mat')];

 avg_SC=zeros(68,68);
for i=1:length(matfiles)  
    name=strcat('/home/shruti/Desktop/Himadri/Shruti/IITG-final-work/Old-Controls/',matfiles(i).name);
    load(name,'SC_cap_agg_bwflav1_norm');
    disp(matfiles(i).name);
    thresh_SC=threshold_proportional(SC_cap_agg_bwflav1_norm,0.11); % 30 percent of strongest weights preserved
    avg_SC=avg_SC+thresh_SC;
    save(name,'thresh_SC','-append');
    clearvars -except avg_SC matfiles i
end
avg_SC_old=avg_SC/length(matfiles);
save('avg_SC_old.mat','avg_SC_old');
 