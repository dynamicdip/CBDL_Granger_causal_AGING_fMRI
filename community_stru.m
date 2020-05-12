clc
clear all
close all
%load('young_old_FC_cc_68.mat','young_old_data');
load('names_data.mat','names_data');
fp = fopen('../controls-gender-age.csv');
out = textscan(fp,'%s%s%s','delimiter',',');
fclose(fp)
% for i=1:49
% FC_ALL = reshape(young_old_data(i,:,:),68,68);
% fisherz = atanh(FC_ALL);
% %    SC_ALL = reshape(young_old_scdata(index,:,:),68,68);
%     thresh_FC_norm= weight_conversion(fisherz,'autofix'); %binarized adjacency matrix
% thresh_FC_norm= weight_conversion(thresh_FC_norm,'normalize');
% thresh_FC=threshold_proportional(thresh_FC_norm,0.1);
%     adjFC= weight_conversion(thresh_FC,'binarize');
% [comm_struct,~] = comm_detect_louvain(thresh_FC_norm,0,0,0,0);
% commStructure{i}=comm_struct;
% commname{i}=names_data{i}.name;
% end

load('result_threshnorm_fc.mat','result');
%for i=1:49
    
for j=1:49
    y(j,1)=max(result{j}.communitymodules);
    y(j,2)=str2double(out{3}(j+1));
end
bar(y(:,1));



