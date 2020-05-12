clear all;
close all;
clc;
%%% YOUNG CONTROLS %%%
matfiles = dir('/home/shruti/Desktop/Himadri/Shruti/IITG-final-work/Young-Controls/*_fMRI_new.mat'); 
avg_FC=zeros(68,68);
adjFC=zeros(68,68);
for i=1:length(matfiles)  
    name=strcat('/home/shruti/Desktop/Himadri/Shruti/IITG-final-work/Young-Controls/',matfiles(i).name);
    load(name,'FC_cc_68');
    disp(matfiles(i).name);
   
    thresh_FC=threshold_proportional(FC_cc_68,0.075);
    adjFC=weight_conversion(thresh_FC,'binarize');
    avg_FC=avg_FC+adjFC;
   save(name,'thresh_FC','-append');
    clearvars -except avg_FC matfiles i adjFC
end
avg_FC_young=avg_FC/length(matfiles);
save('avg_FC_young.mat','avg_FC_young');

%%% OLD CONTROLS %%%%

matfiles = dir('/home/shruti/Desktop/Himadri/Shruti/IITG-final-work/Old-Controls/*_fMRI_new.mat'); 
avg_FC=zeros(68,68);
adjFC=zeros(68,68);
for i=1:length(matfiles)  
    name=strcat('/home/shruti/Desktop/Himadri/Shruti/IITG-final-work/Old-Controls/',matfiles(i).name);
    load(name,'FC_cc_68');
    disp(matfiles(i).name);
   
    thresh_FC=threshold_proportional(FC_cc_68,0.075);
    adjFC=weight_conversion(thresh_FC,'binarize');
    avg_FC=avg_FC+adjFC;
   save(name,'thresh_FC','-append');
    clearvars -except avg_FC matfiles i adjFC
end
avg_FC_old=avg_FC/length(matfiles);
save('avg_FC_old.mat','avg_FC_old');