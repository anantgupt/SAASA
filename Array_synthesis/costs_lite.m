%% This function finds all cost function in 2D data
% peaks: local maxima
% bw: 3dB beamwidth of mainlobe(double sided in U, double sided in V)
% id: Peak location(double sided in U, double sided in V)
% Data should be in dB domain.
%%
function [msll,indx2]= costs_lite(pattern,mask2)
mask=true(5);mask(13)=0;
p=pattern.*(pattern>ordfilt2(pattern,24,mask)).*mask2;
[pm,indx]=max(p(:));
p(indx)=-90;
[pm2,indx2]=max(p(:));
msll=pm2-pm;
end