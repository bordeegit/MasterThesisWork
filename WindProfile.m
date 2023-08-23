function [W] = WindProfile(w_r,phi_r,h_r,surf_rough,H)
% WINDPROFILE Outputs the wind vector in XYZ at height(s) H
%   The output vector has w_x, w_y, w_z stacked on top of each other, hence
%   if the input H is a vector of height, the output will need to be parsed
%   into the single vectors, based on the size of H

% --- SURFACE WIND PROFILE --- %
w_r_x = w_r*cos(phi_r);
w_r_y = w_r*sin(phi_r);

w_x = w_r_x*(log(H/surf_rough)/log(h_r/surf_rough))';
w_y = w_r_y*(log(H/surf_rough)/log(h_r/surf_rough))';
w_z = repelem(0,size(H,2))';

W = [w_x;w_y;w_z];

