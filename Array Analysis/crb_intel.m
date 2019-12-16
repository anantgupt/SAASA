function [CRBi,CRBu,CRBv,ZZBi,ZZBu,ZZBv,curv]=crb_intel(txf,tyf,snrdr,N)
snra=snrdr;colr='rbmycgk';colr=repmat(colr,[1,5]);
mkr={'-','--','-.'};mkr=repmat(mkr,[1,10]);
%% CRB for s(t)=1 model
ifish=(2*pi)^2*[sum((txf-mean(txf)).^2),sum((txf-mean(txf)).*(tyf-mean(tyf)));sum((txf-mean(txf)).*(tyf-mean(tyf))),sum((tyf-mean(tyf)).^2) ];
% ifish=[mean((txf-mean(txf)).^2),mean((txf-mean(txf)).*(tyf-mean(tyf)));mean((txf-mean(txf)).*(tyf-mean(tyf))),mean((tyf-mean(tyf)).^2) ];
Jp=1.3426*eye(2);
% Jp
sigma=(10.^(-snra/10));
% CRB=inv(ifish);
% CRBi=sqrt(eigs(CRB,1)/2)*sqrt(sigma);CRBi=sqrt(sqrt(det(CRB))/2)*sigma;
% CRBu=sqrt(CRB(1,1)/2)*sqrt(sigma);
% CRBv=sqrt(CRB(2,2)/2)*sqrt(sigma);
%% Covariance ellipse
NP = 80;
alpha  = 2*pi/NP*(0:NP);
circle = [cos(alpha);sin(alpha)];
ns = sqrt(6);%From chi square table for p=.95

%% Bayesian CRB, Contribution of Proir J_P
for j=1:length(snra)
    Jf=ifish*2/sigma(j);
    CRB=inv(Jf*N+Jp);

    CRBi(j)=sqrt(det(CRB));
    C=chol(CRB,'lower');
    ellip = ns*C*circle;
%         CRBu(j)=ellip(1,1)/ns;%WARNING sqrt(CRB(1,1));
%     CRBv(j)=ellip(2,NP/4+1)/ns;%WARNINGsqrt(CRB(2,2));
[CRBu(j),ind_max]=max(sqrt(diag(circle'*CRB*circle)));
[CRBv(j),ind_min]=min(sqrt(diag(circle'*CRB*circle)));
%     curv(j)=plot(ellip(1,:),ellip(2,:),['--',colr(j)]);hold on;
%     compass(ellip(1,ind_max),ellip(2,ind_max),['--',colr(j)]);
%         compass(ellip(1,ind_min),ellip(2,ind_min),['--',colr(j)]);
    curv(j)=plot(ellip(1,:),ellip(2,:),['k',mkr{j}]);hold on;
    compass(ellip(1,ind_max),ellip(2,ind_max),['k',mkr{j}]);
        compass(ellip(1,ind_min),ellip(2,ind_min),['k',mkr{j}]);
%       CRBu(j)=sqrt(trace(CRB));CRBv(j)=sqrt(det(CRB));%WARNING: Remove this 
for c=1:NP+1
    t_proj=txf*circle(1,c)-tyf*circle(2,c);
    CRBe(j,c)=sigma(j)/2/((2*pi)^2*N*sum((t_proj-mean(t_proj)).^2));
end
ellipe=[ns*[circle(1,:).*sqrt(CRBe(j,:))];ns*[circle(2,:).*sqrt(CRBe(j,:))]];
    plot(ellipe(1,:),ellipe(2,:),['g',mkr{j}]);
end

%% Based on Van trees paper
D=[txf(:),tyf(:)]';
M=length(txf);
gamma=(M./sigma).^2./((sigma+M)./sigma);
J=N*2/M*(2*pi)^2*D*(eye(M)-ones(M)/M)*D';
CRB=inv(J);
%% Closed form ZZB
Ru=eye(2)/3;
for j=1:length(snra)
ZZB=Ru*2*qfunc(sqrt(N*gamma(j)/2))+CRB./gamma(j).*gammainc(N*gamma(j)/4,1.5);
ZZBi(j)=sqrt(sqrt(det(ZZB)));
ZZBu(j)=sqrt(ZZB(1,1));
ZZBv(j)=sqrt(ZZB(2,2));
end
%% CRB from s(t) is random model
% CRBi=sqrt(eigs(CRB,1)./gamma);CRBi=sqrt(sqrt(det(CRB))./gamma);
% CRBu=sqrt(CRB(1,1)./gamma);
% CRBv=sqrt(CRB(2,2)./gamma);
%% Plotting
% figure(121)
% hold on;
% subplot(1,3,1)
% semilogy(snra,ZZBi,'o-');hold on;
% semilogy(snra,CRBi,'x-');hold on;
% xlabel('SNR(dB)');ylabel('MSE');
% legend('ZZB');
% title('Error Bounds for Range');
% subplot(1,3,2)
% semilogy(snra,ZZBu,'o-');hold on;
% semilogy(snra,CRBu,'x-');hold on;
% xlabel('SNR(dB)');ylabel('MSE');
% legend('ZZB');
% title('Error Bounds for Range');
% subplot(1,3,3)
% semilogy(snra,ZZBv,'o-');hold on;
% semilogy(snra,CRBv,'x-');hold on;
% xlabel('SNR(dB)');ylabel('MSE');
% legend('ZZB');
% title('Error Bounds for Range');
end