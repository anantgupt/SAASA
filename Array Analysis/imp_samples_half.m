function [posa,M_lat,indices]=imp_samples_half(Nt,Rc_u)
% Nt: Desired number of samples
% R_c: sin^2(\theta), theta is elevation from broadside.

fun_u=@(v,u)1./(2*pi*(1-sqrt(1-Rc_u))*sqrt(ones(size(u))-u.^2-v.^2));

if Nt==1
    FU=0;FV=0;Fw=1;
else
    %% Sampling from uniform grid
%     r=linspace(0,sqrt(Rc_u),Nt);theta=linspace(0,1440,Nt);
%     U1=linspace(-sqrt(Rc_u),sqrt(Rc_u),round(sqrt(Nt)));V1=U1;
%     [Ug,Vg]=meshgrid(U1,V1);
%     mask3=find((Ug.^2+Vg.^2) < Rc_u);
%     FU=Ug(mask3);FV=Vg(mask3);Nt=length(mask3);
% % evaluate fun for each (u,v) sample
% target = fun_u(FU, FV);
% calculate importance weight
% Fw = target .*Nt;
    %% Normalizing unit sphere points inside unit cube
%     x=linspace(-sqrt(Rc_u),sqrt(Rc_u),round(sqrt(Nt)));z=x(x>0);
%     [xg,yg,zg]=meshgrid(x,x,z);
%     rg=sqrt(xg.^2+yg.^2+zg.^2);
%     ids=rg<1 & (xg.^2+yg.^2)./(rg.^2)<Rc_u;
%     FU=xg(ids)./rg(ids);FV=yg(ids)./rg(ids);Nt=length(FU);
    %% Using algorithm in Wolfram Mathworld
%      x=linspace(-sqrt(Rc_u),sqrt(Rc_u),round(sqrt(Nt)));
%      [xg,yg]=meshgrid(x,x);
%      rg=xg.^2+yg.^2;
%      ids=rg<1;
%      FU=2*xg(ids).*sqrt(1-rg(ids));FV=2*yg(ids).*sqrt(1-rg(ids));
     %% Using equal area approx
     a=pi*(1-sqrt(1-Rc_u))/Nt;d=sqrt(a);n=1;
     M_lat=round(asin(sqrt(Rc_u))/d);indices=ones(M_lat+1,1);
     d_lat=asin(sqrt(Rc_u))/M_lat;
     d_lon=a/d_lat;
     for i=1:M_lat
        lat_a=(i-0.5)*asin(sqrt(Rc_u))/M_lat;
        M_lon=round(pi*sin(lat_a)/d_lon);
        nn=n+M_lon-1;
        FU(n:nn)=sin(lat_a)*cos(((1:M_lon)-0.5)*pi/M_lon);% Modified to sample only halfcircle
        FV(n:nn)=sin(lat_a)*sin(((1:M_lon)-0.5)*pi/M_lon);
        n=nn+1;
        indices(i+1)=n;
     end
     Fw=ones(size(FU));
end

%% Set Output
Fw = Fw ./ sum(Fw);
posa(1,:) = FU;
posa(2,:) = FV;
posa(3,:)=Fw;
end