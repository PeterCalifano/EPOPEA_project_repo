%% Set Latex font
set(groot,'defaulttextinterpreter','latex')
set(groot,'defaultLegendInterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex')
%%
clearvars; clc; close all

% TCS orbiter Transient

time_FB = load('TimeVectorFB.mat');
Eclipse_index = load('Eclipse_indecesFB.mat');
Angles_FB = load('AnglesFB.mat');

% Retrieve data 2nd flyby
i_zen = Angles_FB.RelPos{2}(:,1); % deg
i_tan = Angles_FB.RelPos{2}(:,2);
i_tran = Angles_FB.RelPos{2}(:,3);
i_out = Angles_FB.RelPos{2}(:,4);
time_2FB = time_FB.time{2}; % second
% Eclipse2 = Eclipse_index{2}; % ---> no eclipses in 2nd fb

% time for transient
t0 = time_2FB(1);
tf = time_2FB(end);

% interpolate
i_zen_tr = @(t) interp1(time_2FB, i_zen, t);


interp1(time_2FB, i_zen, 10000)