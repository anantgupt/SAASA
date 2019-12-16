%% This function finds all cost function in 2D data
% peaks: local maxima
% bw: 3dB beamwidth of mainlobe(double sided in U, double sided in V)
% id: Peak location(double sided in U, double sided in V)
% Data should be in dB domain.
%%
function [cost,msll,beamw,eccen,GD]= costs_lite2(pattern,mask2)
global SU;global SV;global pbi;
alpha=0.2;beta=0.5;gamma=0.5;psi=1;%20
mask=true(3);mask(5)=0;
p=pattern.*(pattern>ordfilt2(pattern,8,mask)).*mask2;
[pm,indx]=max(p(:));
p(indx)=-90;
[pm2,~]=max(p(:));
% scale by square root aperture area
%     ang=linspace(0,pi*2*(1-1/Nres),Nres);
%     txf=[tx+0.75,tx+0.75,tx-0.75,tx-0.75];tyf=[ty+0.9,ty-0.9,ty+0.9,ty-0.9];
%     ap=zeros(1,Nres);
%     for i=1:Nres
%         apt=(txf+1j*tyf)*exp(-1j*ang(i));
%     ap(i)=max(real(apt))-min(real(apt));
%     end
%         area_ap=trapz((2./ap).^2)*(ang(2)-ang(1))/2;% Area of inverse aperture in UV
%% Correct Beamwidth for fan beams based on area of ellipse, 
indxu=find(pattern>(pm-3));indxd=find(pattern<(pm-3) & pattern>(pm-5));
theta=asin(sqrt(SU(indx)^2+SV(indx)^2));
% [amaj,~]=max((SU(indxu)-SU(indx)).^2+(SV(indxu)-SV(indx)).^2);
% [amin,~]=min((SU(indxd)-SU(indx)).^2+(SV(indxd)-SV(indx)).^2);
% bd1=180/pi*asin(sqrt(amaj)/cos(theta))*2;
% bd2=180/pi*asin(sqrt(amin)/cos(theta))*2;
%% Aperture inverse beamwidth
% area_bp=pi*sqrt(amaj*amin)/(cos(theta)^2);% Area of 3dB beam in UV
% scale_fac=sqrt(area_ap/area_bp); % Scale factor, see Feb 15 ppt
%     fill(20./ap.*cos(ang),20./ap.*sin(ang),'w','linewidth',2);axis([-10,10,-10,10]);grid on;hold on;%Inverse Aperture
% scatter(10*scale_fac*SU(indxu),10*scale_fac*SV(indxu),[],pattern(indxu),'*');view(2);colorbar;%3dB beam scaled
%msll=pm2-pm-beta*log10(area)+gamma*(max(ap)-min(ap));
% beamw=area_ap*scale_fac^2;
%% ALternative beamwidth computation using 3D dot product
[amaj,~]=min(SU(indxu)*SU(indx)+SV(indxu)*SV(indx)+cos(theta)*sqrt(1-SU(indxu).^2-SV(indxu).^2));%acos is decreasing function
[amin,~]=max(SU(indxd)*SU(indx)+SV(indxd)*SV(indx)+cos(theta)*sqrt(1-SU(indxd).^2-SV(indxd).^2));%acos is decreasing function
bd1=acos(amaj)*2*180/pi;
bd2=acos(amin)*2*180/pi;
beamw=bd1;% was bd1*bd2; May16, 2018
msll=(pm2-pm);
% eccen=sqrt(max(ap)^2-min(ap)^2)/max(ap);
eccen=sqrt(bd1^2-bd2^2)/bd1;
GD=gd(db2mag(pattern+pbi),SU,SV);
cost=ov_cost(beamw,msll,eccen,GD);% eccen=scale_fac;
end