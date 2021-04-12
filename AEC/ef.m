function [Properties,Input,Others]=ef(varargin)

clc
efs=fopen('sampling.txt');
sampling_scan=textscan(efs,'%s %f','delimiter','=');
fclose(efs);
sampling_variables={'U0','Q','d0','SOT','theta_s','theta_v'};
for i_check=1:length(sampling_variables)
ind(i_check)=isempty(find(strcmpi(varargin,sampling_variables{i_check})==1,1));
end

ind_update=find(ind==0);%the input variable

for i_up=1:length(sampling_variables)
sampling_inputs(i_up)=sampling_scan{2}(strcmpi(sampling_scan{1},sampling_variables{i_up}));
end

for i_up=1:length(ind_update)
ind_var(i_up)=  find(strcmpi(varargin,sampling_variables{ind_update(i_up)})==1); % find the position of the variable that's to be updated from varargin
sampling_inputs(ind_update(i_up))=varargin{ind_var(i_up)+1};
end
end