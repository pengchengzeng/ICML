clear
clc
close all

%load data. Take dataset1 for example
load('./dataset1/D11.mat');
load('./dataset1/D21.mat');
load('./dataset1/D12.mat');
load('./dataset1/D22.mat');
load('./dataset1/labels.mat'); 

%set the default value for Kij and N; we can also tune these
%hyperperameters based on the method by coupleCoC+.
Kij=10; N = 9; alpha = 1; nh = Kij; iter = 5;
nrowcluster11=Kij;nrowcluster12=Kij;nrowcluster21=Kij;nrowcluster22=Kij;
ncolcluster=N;

%Initialization. 
%Method 1: random initialization as follows:
%rng(1);
%Cx110 = randsample(nrowcluster11, size(p11,1),true); 
%Cx120 = randsample(nrowcluster12, size(p12,1),true); 
%Cx210 = randsample(nrowcluster21, size(p21,1),true); 
%Cx220 = randsample(nrowcluster22, size(p22,1),true); 
%Cy0 = randsample(ncolcluster, size(p11,2),true); 
%Method 2: Get the initialization based on ITCC for each view of data. 
%Here we directly give the initialized values.
load('./dataset1/Cx110.mat');
load('./dataset1/Cx120.mat');
load('./dataset1/Cx210.mat');
load('./dataset1/Cx220.mat');
load('./dataset1/Cy0.mat');


%scICML algorithm
[Cy, ~, ~, ~, ~, ~] = scICML_main(p11,p12,p21,p22,Cy0, Cx110, Cx120, Cx210, Cx220, alpha,iter,nh);

%present the result
NMIARI(Cy, gnd);

