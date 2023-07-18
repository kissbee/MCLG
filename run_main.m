
clc;
clear memory;
clear;
clear all;
warning('off')
addpath('data')
addpath('func')

% make dir to save results 
file_path_name_ = 'result';
if exist(file_path_name_,'file')==0   %if not exist,make one
    mkdir(file_path_name_);
end

dataname = 'BRCA';                  % select dataset
disp(dataname);                     
load(dataname);                     % load data

Test_MCLG(dataname,X,Y);            % test data dimension: n * dv

