% For numerical, uniform in UV
function [ZZBi,ZZBu,ZZBv,curv]=zzb_intel7(snrdr,txf,tyf,N,ha1,vva1,A1,Ns,pb)
% txf,tyf must be single configuration 128x1
snra=snrdr;M=length(txf);Rc_delta=1;Nd=512;NP=60;Rc_u=0.25;
%% For ZZB Pe Numerical computation
fun_delta=@(v,u)1./(2*pi*(1-sqrt(1-Rc_delta))*sqrt(ones(size(u))-u.^2-v.^2));
colr='rbmycgk';colr=repmat(colr,[1,5]);
mkr={'-','--','-.'};mkr=repmat(mkr,[1,10]);
%% Do importance sampling for only broadside (posa bfore)
posa=imp_samples(Ns,Rc_u);% was [0;0;1]
Nres=512;
U = linspace(0,1-1/Nres,Nres).^2;V = linspace(-1,1,Nd+1);%.*abs(linspace(-1,1,Nd+1));
[SU,SV] = meshgrid(U,V);
mask2=(SU.^2+SV.^2) < 1;
indx=find((SU.^2+SV.^2)>1);
%% Numerical computation of Pe
ha=U;
alpha  = pi/NP*(0:NP-1);c95=sqrt(6);
circle = [cos(alpha);sin(alpha)];
circle2 = -circle;
% figure(30)
% hold off;
for c=1:NP
    VansU=zeros(length(ha),1);
    %     SU2=SU*circle(1,c)-SV*circle(2,c);SV2=SU*circle(2,c)+SV*circle(1,c);
    %     pattern=bpmra_lite(txf-mean(txf),tyf-mean(tyf),SU2,SV2,zeros(size(txf)),indx);
    %     pattern2=bpmra_lite(txf,tyf,SU3,SV3,zeros(size(txf)),indx);
    tx_proj=txf*circle(1,c)-tyf*circle(2,c);
    ty_proj=txf*circle(2,c)+tyf*circle(1,c);
    %     pattern=bpmra_lite(tx_proj,ty_proj,SU,SV,zeros(size(txf)),indx);
    pattern=bpmra_lite(tx_proj-mean(tx_proj),ty_proj-mean(ty_proj),SU,SV,zeros(size(txf)),indx);
    %     figure(29)
    % %     subplot(2,1,1)
    %     hold off;
    %     surf(SU,SV,max(real(pattern),0),'linestyle','none');view(2);colorbar;hold on;
    % %     surf(SU3,SV3,max(real(pattern2),0),'linestyle','none');view(2);
    % %     quiver(0,0,circle(1,c),circle(2,c));
    % %     subplot(2,1,2)
    % %     surf(SU3,SV3,max(real(pattern2),0),'linestyle','none');view(2);colorbar;hold on;
    if 0
        for h=1:length(ha)
            % %         subplot(2,1,1)
            %         [mval,mid]=max(real(pattern(:,h)));
            %         plot3(SU(mid,h),SV(mid,h),mval,'r.');
            for i=1:length(snra)
                sigma=(10.^(-snra(i)/10));
                %             Q(i,h)=qfunc(sqrt((M-max(real(pattern(:,h))))*N/sigma));
                % Non-coherent case
                if 0 % costly
                    a=sqrt((1-sqrt(1-abs(pattern(:,h)/M).^2))*M*N/4/sigma);
                    b=sqrt((1+sqrt(1-abs(pattern(:,h)/M).^2))*M*N/4/sigma);
                    Q(i,h)=max(marcumq(a,b))-0.5*exp(-N*M/4/sigma)*besseli(0,max(abs(pattern(:,h)))*N/4/sigma);
                    
                    
                else %fast-using fact that Q0, I0 are monotone
                    a=sqrt((1-sqrt(1-max(abs(pattern(:,h)/M)).^2))*M*N/2/sigma);
                    b=sqrt((1+sqrt(1-max(abs(pattern(:,h)/M)).^2))*M*N/2/sigma);
                    Q(i,h)=marcumq(a,b)-0.5*exp(-(a^2+b^2)/2)*besseli(0,a*b);
                end
            end
        end
        %     ZZBu1(c)=trapz(ha,(VansU.').*ha);
        for i=1:length(snra)
            ZZBu1(i,c)=trapz(ha,valley(Q(i,:).').*ha);
            %         plot(valley(Q(i,:).'));hold on;
        end
        %     hold off;
    else % Analytical ZZB: 
                V=linspace(-0.5,0.5,65);%.*abs(linspace(-1,1,129));
        pat_fun=@(h)max(abs(bpmra_lite(tx_proj-mean(tx_proj),ty_proj-mean(ty_proj),ones(size(V)).'*h,...
            V.'*ones(size(h)),zeros(size(txf)),(ones(size(V)).'*h).^2+(V.'*ones(size(h))).^2>1)));
                    

        for i=1:length(snra)
            sg=(10.^(-snra(i)/10));
            af=@(h)sqrt((1-sqrt(1-abs(pat_fun(h)/M).^2))*M*N/2/sg);
            bf=@(h)sqrt((1+sqrt(1-abs(pat_fun(h)/M).^2))*M*N/2/sg);
            Qfun2=@(h)h.*valley(marcumq(af(h),bf(h),1)-0.5*exp(-(af(h).^2+bf(h).^2)/2).*besseli(0,af(h).*bf(h)));%abs(pat_fun(h))/4/sigmaj
            ZZBu1(i,c)=integral(Qfun2,0,1,'AbsTol',1e-5);
        end
    end
end
for i=1:length(snra)
    % Plot directly
    ellip=[c95*[circle(1,:).*sqrt(ZZBu1(i,:)),circle2(1,:).*sqrt(ZZBu1(i,:)),circle(1,1).*sqrt(ZZBu1(i,1))]...
        ;c95*[circle(2,:).*sqrt(ZZBu1(i,:)),circle2(2,:).*sqrt(ZZBu1(i,:)),circle(2,1).*sqrt(ZZBu1(i,1))]];
    [ZZBu(i),ind_max]=max(sqrt(ZZBu1(i,:)));
    [ZZBv(i),ind_min]=min(sqrt(ZZBu1(i,:)));
    % Compute R
    %     A=[circle(1,:).^2;2*circle(1,:).*circle(2,:);circle(2,:).^2];
    %     Rel=ZZBu1(i,:)/A;
    %     Rzzb=[Rel(1),Rel(2);Rel(2),Rel(3)];
    %     C=chol(Rzzb,'lower');
    %     ellip = c95*C*circle;
    %     [ZZBu(i),ind_max]=max(sqrt(diag(circle'*Rzzb*circle)));
    %     [ZZBv(i),ind_min]=min(sqrt(diag(circle'*Rzzb*circle)));
    
    curv(i)=plot(ellip(1,:),ellip(2,:),['b',mkr{i}]);hold on;
    compass(ellip(1,ind_max),ellip(2,ind_max),['b',mkr{i}]);
    compass(ellip(1,ind_min),ellip(2,ind_min),['b',mkr{i}]);
    
    %     ZZBv(i)=trapz(ha,valley(VansV.').*ha);
    % A(h)P(h) curves
    %     loglog(ha,(VansU),':x','Color',colr(i));hold on;
    %     loglog(ha,(VansV),'-d','Color',colr(i));hold on;
    % if i==4
    % %     surf(hm,vm,PU,'linestyle','none');
    %     plot(PV.');
    % end
end
% Plotting 3dB BW
% pm=max(real(pattern(:)));
% indxu=find(pattern>(pm/sqrt(2)));
% scatter(SU2(indxu),SV2(indxu),'g.');hold on;
% indxu=find(pattern2>(pm/sqrt(2)));
% scatter(SU3(indxu),SV3(indxu),'m.');
ZZBi=sqrt(ZZBu.*ZZBv);
end