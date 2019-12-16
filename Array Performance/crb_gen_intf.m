function [CRBi,CRBu,CRBv,CSCRBi,CSCRBu,CSCRBv]=crb_gen_intf(txf,tyf,snrdr,N,Phi,g,u,v)
snra=snrdr;colr='rbmycgk';colr=repmat(colr,[1,5]);
mkr={'-','--','-.'};mkr=repmat(mkr,[1,10]);
txf=txf-mean(txf);
tyf=tyf-mean(tyf);
%% CRB for s(t)=1 model
part_der=zeros(length(txf),length(g));
for id=1:length(g)
    part_der(:,2*id-1)=g(id)*(2*pi*txf.*exp(1j*2*pi*(txf*u(id)+tyf*v(id))));
    part_der(:,2*id)=g(id)*(2*pi*tyf.*exp(1j*2*pi*(txf*u(id)+tyf*v(id))));
end
ifish=real(part_der'*part_der);
%% Add Interference (only for 2 Target)
% ifish(1:2,3:4)=0.5*ifish(1:2,3:4);% Half off diagonal 2x2 block
% ifish(3:4,1:2)=0.5*ifish(3:4,1:2);
% ifish(1:2,1:2)=ifish(1:2,1:2)+ifish(1:2,3:4);% Add off diagonal 2x2 block to diagonal
% ifish(3:4,3:4)=ifish(3:4,3:4)+ifish(3:4,1:2);

part_derCS=Phi*part_der;
ifishCS=real(part_derCS'*part_derCS);
%% Add interference
% ifishCS(1:2,3:4)=0.5*ifishCS(1:2,3:4);% Half off diagonal 2x2 block
% ifishCS(3:4,1:2)=0.5*ifishCS(3:4,1:2);
% ifishCS(1:2,1:2)=ifishCS(1:2,1:2)+ifishCS(1:2,3:4);% Add off diagonal 2x2 block to diagonal
% ifishCS(3:4,3:4)=ifishCS(3:4,3:4)+ifishCS(3:4,1:2);

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
%% Bayesian CRB, Contribution of Proir J_P
for j=1:length(snra)
    Jf=ifish*2/sigma(j);
    CRB=inv(Jf*N);%+Jp);
    
    Jf2=ifishCS*2/sigma(j)*(size(Phi,1)/size(Phi,2));%Noise scales bcoz Phi is normalized to preserve signal magnitude
    %     Jf2=ifishCS*2/sigma(j)*(size(Phi,1)/size(Phi,2));
    
    CRB2=inv(Jf2*N);
    %%
    %     C=chol(CRB,'lower');
    %     ellip = ns*C*circle;
    % %         CRBu(j)=ellip(1,1)/ns;%WARNING sqrt(CRB(1,1));
    % %     CRBv(j)=ellip(2,NP/4+1)/ns;%WARNINGsqrt(CRB(2,2));
    % [CRBu(j),ind_max]=max(sqrt(sum(ellip.^2,1))/ns);
    % [CRBv(j),ind_min]=min(sqrt(sum(ellip.^2,1))/ns);
    %     curv(j)=plot(ellip(1,:),ellip(2,:),['k',mkr{j}]);hold on;
    %     compass(ellip(1,ind_max),ellip(2,ind_max),['k',mkr{j}]);
    %         compass(ellip(1,ind_min),ellip(2,ind_min),['k',mkr{j}]);
    %%
    for id=1:length(g)
        CRBl1=CRB(2*id-1:2*id,2*id-1:2*id);
        CRBi(j,id)=sqrt(det(CRBl1));
        [CRBu(j,id),~]=max(sqrt(diag(circle'*CRBl1*circle)));
        [CRBv(j,id),~]=min(sqrt(diag(circle'*CRBl1*circle)));
        %         CRBu(j,id)=sqrt(trace(CRBl1));CRBv(j,id)=sqrt(det(CRBl1));%WARNING: Remove this
        CRBl2=CRB2(2*id-1:2*id,2*id-1:2*id);
        CSCRBi(j,id)=sqrt(det(CRBl2));
        %         CSCRBu(j,id)=sqrt(trace(CRBl2));CSCRBv(j,id)=sqrt(det(CRBl2));
        [CSCRBu(j,id),~]=max(sqrt(diag(circle'*CRBl2*circle)));
        [CSCRBv(j,id),~]=min(sqrt(diag(circle'*CRBl2*circle)));
    end
end
 CSCRBu=sqrt(CSCRBu.^2+CSCRBv.^2);%WARNING: u is overall error (MSE)
    CRBu=sqrt(CRBu.^2+CRBv.^2);%WARNING: u is overall error (MSE)
end