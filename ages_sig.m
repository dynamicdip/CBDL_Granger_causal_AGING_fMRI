clc
clear all
close all
load('granger_allsubjects_sig.mat','subject_sig');
load('names_data.mat','names_data');
fp = fopen('../controls-gender-age.csv');
out = textscan(fp,'%s%s%s','delimiter',',');
fclose(fp)
age_sigsum(1:3,:,:)=zeros(3,10,10);
age_data_count(1:3,:)=zeros(3,1);
for i=1:49
    age = str2double(out{3}(i+1));
    if age <= 25
        age_sigsum(1,:,:)=age_sigsum(1,:,:)+subject_sig(i,:,:);
        age_data_count(1,:) = 1 + age_data_count(1,:);
    elseif age > 25 &&  age <= 50
        age_sigsum(2,:,:)=age_sigsum(2,:,:)+subject_sig(i,:,:);
        age_data_count(2,:) = 1 + age_data_count(2,:);
    else
        age_sigsum(3,:,:)=age_sigsum(3,:,:)+subject_sig(i,:,:);
        age_data_count(3,:) = 1 + age_data_count(3,:);
    end
end

for i=1:3
    age_sigsum_result(i,:,:) = age_sigsum(i,:,:)./age_data_count(i);
end

