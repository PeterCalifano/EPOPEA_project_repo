function [alpha,alpha_dot]=pointingcam(h,FOV,Np,texp,dx)
%% PROTOTYPE
% [alpha,alpha_dot]=pointingcalc(h,FOV,Np,texp,dx)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Computes the required pointing budget of an optical instrument
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% h [km] - altitude 
% FOV [°] - FOV of the instrument
% Np [-] - number of pixels per axis of the instrument
% texp [s] - exposure time
% dx [m] - spatial resolution
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% alpha [°] - pointing accuracy requirement
% alpha_dot [°/s] - pointing stability requirement
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 21/11/22 - Fabrizio Maccari - first implementation
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades


l=2*h*tan(deg2rad(FOV/2))*1e+3; %[m] - swath width
p=l/Np; %[m] - on ground size of a pixel

Dl=l/3; %accepted drift of the centre of the image

alpha=rad2deg(atan(Dl/(h*1e+3))); %[°] pointing accuracy

p_drift=floor(dx/p); %accepted number of dislocated pixels
res=p_drift*p; %resolution
alpha_d=rad2deg(atan(res/(h*1e+3))); %drift
alpha_dot=alpha_d/texp;



