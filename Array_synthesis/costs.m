%% This function finds all cost function in 2D data
% peaks: local maxima
% bw: 3dB beamwidth of mainlobe(double sided in U, double sided in V)
% id: Peak location(double sided in U, double sided in V)
% Data should be in dB domain.
%%
function [msll,asll,bwa,msls,ptge,umax,vmax,id,id2,gd2,pm,pm2,ecc]= costs(pattern,U,V,u0,v0,GD)
mask=true(3);mask(5)=0;N=size(pattern,2);midi=N/2+1;
mask2=(U.^2+V.^2) < 1-2/N;% Removes unnecessary peaks found on boundary
p=pattern.*(pattern>ordfilt2(pattern,8,mask)).*mask2;
[pm,indx]=max(p(:));
id(2)=ceil(indx/N);
id(1)=indx-(id(2)-1)*N;
p(isnan(p))=0;% NaN were encountered for + config bcoz -Inf.
[peaks,id2t]=sort(p(:),'descend');
indx2=id2t(2);
id2(2)=ceil(indx2/N);
id2(1)=indx2-(id2(2)-1)*N;
umax=U(indx);vmax=V(indx);
ptge=sqrt(((umax-u0)^2+(vmax-v0)^2));
msls=(U(indx2)-U(indx))^2+(V(indx2)-V(indx))^2;
theta=asin(sqrt(U(indx)^2+V(indx)^2));
% bw(1)=180/pi*asin((find(pattern(N/2+1:end,id(2))<(pm-3),1)/cos(theta))*2/N)*2;%mainbeam
% bw(2)=180/pi*asin((find(pattern(id(1),N/2+1:end)<(pm-3),1)/cos(theta))*2/N)*2;
% bw(1)=180/pi*asin((find(pattern(N/2+1:end,midi)<(pm-3),1)/cos(theta))*2/N)*2;%broadside
% bw(2)=180/pi*asin((find(pattern(midi,N/2+1:end)<(pm-3),1)/cos(theta))*2/N)*2;
peaks=peaks(1:find(peaks==0,1));
% Interference costs
%% Final Output
msll=peaks(2)-peaks(1);
pm2=peaks(2);
gd2=pm-pattern(round(v0*N/2)+midi,round(u0*N/2)+midi)+GD;
asll=0.5*mag2db(mean(db2mag(peaks*2)))-peaks(1);
% bwa=prod(bw);
%         sllhist=histcounts(peaks-peaks(1),-45:0);
% Beamwidth for fan beams based on area of ellipse in UV
indxu=find(pattern>(pm-3));indxd=find(pattern<(pm-3) & pattern>(pm-5));
% [amaj,imaj]=max((U(indxu)-umax).^2+(V(indxu)-vmax).^2);
% [amin,imin]=min((U(indxd)-umax).^2+(V(indxd)-vmax).^2);
% bd1=180/pi*asin(sqrt(amaj)/cos(theta))*2;
% bd2=180/pi*asin(sqrt(amin)/cos(theta))*2;
%% Alternative beamwidth computation based on 3D dot product
[amaj,~]=min(U(indxu)*U(indx)+V(indxu)*V(indx)+cos(theta)*sqrt(1-U(indxu).^2-V(indxu).^2));%acos is decreasing function
[amin,~]=max(U(indxd)*U(indx)+V(indxd)*V(indx)+cos(theta)*sqrt(1-U(indxd).^2-V(indxd).^2));%acos is decreasing function
bd1=acos(amaj)*2*180/pi;
bd2=acos(amin)*2*180/pi;
%%
bwa=sqrt(bd1.^2+bd2.^2);% NOTE: This is overall BW was bd1*bd2;
ecc=sqrt(bd1^2-bd2^2)/bd1;
% sllhist=[U(indxu(imaj)),V(indxu(imaj)),U(indxd(imin)),V(indxd(imin))];
% figure
%         p=pattern.*(pattern<ordfilt2(pattern,3,mask)).*mask2;% Peak to null beamwidth
% imshow(p./pm);
%% Prob of outlier DML case(Athley2005)
% sf=
% pdml=marcumq(sqrt(sf/2*(1-sqrt(1-rf.^2))),sqrt(sf/2*(1+sqrt(1-rf.^2))))...
%     +exp(-sf/2)*(besselj(rf*sf/2)-2^(-2*Nf+1).*(besselj(rf*Sf/2)*sum(nchoosek(2*Nf-1,0:(Nf-1)))-sum(besseli(1:(Nf-1),rf*Sf/2)).*().*sum(nchoosek(2*Nf-1,0:(Nf-1-m))) ));
end