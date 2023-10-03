function [W] = SmoothWindProfile(p_x,p_y,p_z, WxFnc, WyFnc)
%SMOOTHWINDPROFILE Produces wind given 3d point (p_x,p_y,p_z) 
%   The wind depends on a cubic function on the z and a asymmetric function
%   wrt x,y 
W = zeros(3,1);

W(1) = feval(WxFnc,p_z);
W(2) = feval(WyFnc,p_x,p_y);

end

