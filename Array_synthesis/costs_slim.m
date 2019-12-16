%% This function finds all cost function in 2D data
% peaks: local maxima
% bw: 3dB beamwidth of mainlobe(double sided in U, double sided in V)
% id: Peak location(double sided in U, double sided in V)
% Data should be in dB domain.
%%
function msll= costs_slim(pattern,mask2)
mask=true(3);mask(5)=0;
p=pattern.*(pattern>ordfilt2(pattern,8,mask)).*mask2;
[pm,indx]=max(p(:));
p(indx)=-90;
[pm2,~]=max(p(:));
msll=pm2-pm;
end