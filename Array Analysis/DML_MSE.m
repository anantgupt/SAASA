function [Ru,Rv,Ri,c]=DML_MSE(snrdr,txf,tyf,N,Ns,stdmm,pb)
Nu=100;M=length(txf);Rc_u=(1/Nu).^2;colr='rgbmkc';colr=repmat(colr,[1,5]);mkr={'-','--','-.'};mkr=repmat(mkr,[1,10]);
% U = [-1:0.85/Nu:-0.02,-exp(-4:-1:-Nu/2),0,exp(-Nu/2:-4),0.02:0.85/Nu:1];
% U=[(-1:1/Nu:-0.17),-(-0.4:1/Nu:0).^2,(1/Nu:1/Nu:0.4).^2,0.17:1/Nu:1];
% U=[(-1:2/Nu:-0.27),-(-0.5:2/Nu:0).^2,(2/Nu:2/Nu:0.5).^2,0.27:2/Nu:1];
U=[-1:2/Nu:-0.51,(-0.5:1/Nu:0.5),0.51:2/Nu:1];
% U=linspace(-0.707,0.707,Nu);
V = U;
[SU,SV] = meshgrid(V,U);
mask2=find((SU.^2+SV.^2) < 1-2/Nu);
txf=txf-mean(txf);tyf=tyf-mean(tyf);
Dict=exp(1j*(2*pi*(txf*(SU(mask2)).'+tyf*(SV(mask2)).')));
NDict=numel(mask2);
Uf=linspace(-1/Nu,1/Nu,20);
% Uf=linspace(-0.08,0.08,50);Uf=Uf.*abs(Uf);
Vf=Uf;
[SUf,SVf] = meshgrid(Vf,Uf);
Dictf=exp(1j*(2*pi*(txf*(SUf(:)).'+tyf*(SVf(:)).')));
NDictf=numel(SUf);
Navg=400; % Number of averaging

%%
NP = 60;
    alpa  = pi/NP*(0:NP);
    circle = [cos(alpa);sin(alpa)];
    ns = sqrt(6);%From chi square table for p=.95
    MLEe=zeros(length(snrdr),NP+1);
%% Noise to element positions
% stdmm=2;%Standard Deviation

a=zeros(M,1);
% Dict=exp(1j*(2*pi*(txf*(SU(mask2)).'+tyf*(SV(mask2)).'))).*(repmat(patchbeam(SU(mask2),SV(mask2)),1,M)).';
Ru=zeros(1,length(snrdr));Rv=Ru;cov=zeros(2,2,length(snrdr));
posa=imp_samples_half(Ns,Rc_u);Ns=size(posa,2);
rng(5);% same randomization across configs
alpha=exp(2*pi*1j*rand(N,Ns,Navg));
rand_noise=(randn(N,M,Ns,Navg)+1j*randn(N,M,Ns,Navg));
% figure(77)
for i=1:length(snrdr)
snr=snrdr(i);
sigma=db2mag(-snrdr(i)).^2;
    %%
    for t=1:Ns
        u0=posa(1,t);
        v0=posa(2,t);
        w=posa(3,t);
        %         uv1=(u0)^2+(v0)^2;Pue=0;Pve=0;
        %         if uv1<Rc_u
        for nave=1:Navg % alpha phase const over averaging
            txf2=txf;%+randn(M,1)*stdmm/5/sqrt(2);
            tyf2=tyf;%+randn(M,1)*stdmm/5/sqrt(2);
            if pb==1
                Xn=patchbeam(u0,v0)*ones(N,1)*(exp(1j*2*pi*(txf2*u0+tyf2*v0)).')+wgn(N,M,-snr,'complex');%+(randn(N,M)+1j*randn(N,M))*sqrt(sigma/2);
            else
                Xn=squeeze(alpha(:,t,nave))*(exp(1j*2*pi*(txf2*u0+tyf2*v0)).')+squeeze(rand_noise(:,:,t,nave))*sqrt(sigma/2);
%                 Xn=ones(N,1)*(exp(1j*2*pi*(txf2*u0+tyf2*v0)).')+(randn(N,M)+1j*randn(N,M))*sqrt(sigma/2);
            end
%             Rxx=conj(Xn'*Xn)/N;
%             for p=1:NDict
%                 a=Dict(:,p);
%                 proj_res(p)=real(a'*Rxx*a);
%             end
%             [~,id1]=min(trace(Rxx)-proj_res/M);
            proj_res=sum(abs(Dict'*Xn.').^2,2);
            [~,id1]=max(proj_res);
            if 1 % Refine or not
                Dictf_shift=exp(1j*(2*pi*(SU(mask2(id1))*txf+SV(mask2(id1))*tyf)));
%                             for p=1:NDictf
%                                 a=Dictf(:,p).*Dictf_shift;
%                                 proj_resf(p)=real(a'*Rxx*a);
%                             end
%                             [~,id2]=min(trace(Rxx)-proj_resf/M);
                proj_resf=sum(abs(Dictf'*diag(Dictf_shift')*Xn.').^2,2);
                [~,id2]=max(proj_resf);
                eu=(SU(mask2(id1))+SUf(id2)-u0);
                ev=(SV(mask2(id1))+SVf(id2)-v0);
            else
                eu=(SU(mask2(id1))-u0);
                ev=(SV(mask2(id1))-v0);
            end
            Ru(i)=Ru(i)+w*sum(eu.^2)/N/Navg;
            Rv(i)=Rv(i)+w*sum(ev.^2)/N/Navg;
            cov(:,:,i)=cov(:,:,i)+w*([eu,ev].'*[eu,ev])/Navg;
            MLEe(i,:)=MLEe(i,:)+w*abs([eu,ev]*circle).^2/Navg;% MLE ellipse
        end
    end
        %% Covariance ellipse
    
    Re=squeeze(cov(:,:,i));
    C=chol(Re,'lower');
    ellip = ns*C*circle;
    
    [Ru(i),ind_max]=max(sqrt(diag(circle'*Re*circle)));
    [Rv(i),ind_min]=min(sqrt(diag(circle'*Re*circle)));
    
    c(i)=plot(ellip(1,:),ellip(2,:),['r',mkr{i}]);hold on;
    compass(ellip(1,ind_max),ellip(2,ind_max),['r',mkr{i}]);
    compass(ellip(1,ind_min),ellip(2,ind_min),['r',mkr{i}]);
    %% MLE ellipse
    ellipe=[ns*[circle(1,:).*sqrt(MLEe(i,:))];ns*[circle(2,:).*sqrt(MLEe(i,:))]];
    plot(ellipe(1,:),ellipe(2,:),['m',mkr{i}]);hold on;
%     Ru(i)=sqrt(trace(Re));
%     Rv(i)=sqrt(det(Re));% WARNING:Remove this
    Ri(i)=sqrt(det(Re));
end
% Ru=sqrt(Ru);Rv=sqrt(Rv);
end