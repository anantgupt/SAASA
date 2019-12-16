% function o_s = Compressive_Estimationf(rade,per_sub,CS,g_array,N_DoA,K_DoA,M_Phi,n_far,snrdr,N_Phi_Mat,Navg,Nref,date_str)
clear all;close all;
set(0,'DefaultFigureWindowStyle','docked');
%%%%%%%%%%%%%%%%  Simulation Parameters  %%%%%%%%%%%%%%%%%%
rade=0; % Rademacher (QPSK), otherwise GRV
per_sub=0; %Subarray compressive or total
CS=0+per_sub; %Compressive
vec =  @(in) in(:);
cfg_sel=2;% Select configuration

K_DoA=2;
snrdr=[-9,0];
% snrdr=[-20:5:-10,-9:2:-1,0:5:5];
Nref=4;
M_Phi = 32; % Number of Compressive measurements
N_Phi_Mat=1;Navg=8;
date_str=['3jun25_',num2str(cfg_sel)];

calcium0={'r-','b-','g-','k-'};calcium1={'r--','b--','g--','k--'};calcium2={'r-s','b-s','g-s','k-s'};
calcium3={'r-.v','b-.v','g-.v','k-.v'};calcium4={'r:x','b:x','g:x','k:x'};calcium5={'r.-','b.-','g.-','k.-'};
N=1;%# time snapshots was 20
Ns=128;% Number of array elements
pb_case=0;
Nsnr=length(snrdr);
%%%%%%%%%%%%%%%%  Scene Parameters  %%%%%%%%%%%%%%%%%%

N_DoA=4096;%was 256
uniform_DoA=0;
N_DoA2=N_DoA/Navg;
N_rings=floor(sqrt(N_DoA2));

if uniform_DoA
    [DoAs,N_rings,ring_start_ind]=imp_samples_half(N_DoA2,0.25);
    N_DoA2=size(DoAs,2);
    N_DoA=N_DoA2*Navg;
    doa_sep=sqrt(sum(DoAs(1:2,ring_start_ind(1:end-1)).^2));
else
    NP=(N_DoA2/N_rings);
    alpha  = pi/NP*(0:NP-1);
    circle = [cos(alpha);sin(alpha);alpha];
    doa_sep=linspace(0.1,0.7,N_rings).^2;
    DoAs=kron(doa_sep,circle);%imp_samples(N_DoA2,0.25);
    N_DoA2=size(DoAs,2);
    N_DoA=N_DoA2*Navg;
    DoAs(3,:)=ones(1,N_DoA2)/N_DoA2;
end
n_far=round(N_DoA2/2);% Index of first element far away
if K_DoA==1
    g_array=(ones(N_DoA,1)).*exp(2*pi*1i*rand(N_DoA,K_DoA));
else
    g_array=(ones(N_DoA,1)*[sqrt(0.25),1]).*exp(2*pi*1i*rand(N_DoA,K_DoA));
end
%%
Np=4;Nsub=8;
img = imread('rfem.png');             %# Load rfem image
norm_compare=@(in)(in-0.5*(max(in)+min(in)))/(max(in)-min(in));
configs={'Tri','crc','lin','B'};
% config={'Triangle','Circle','Linear','Benchmark'};
config={'BW,MSLL','Benchmark','BW only','All'};%,'Linear','Benchmark'};
Alg_name={'NOMP A','NOMP B','Lasso','NLasso'};
% Alg_name={'NOMP','BPDN (SPGL1)','LASSO (CVX)'};
colr=['r','b','g','k'];
if rade, batch_name='QPSK'; else batch_name='gaus';end
if per_sub==1 dir_name=[date_str,'_CS',num2str(CS),'_persubarray/']; else dir_name=[date_str,'_CS',num2str(CS),'/']; end
eps =0.001;rho=1.0;
%% Subarray beam pattern
Nxr=4;Nyr=4;dolx=0.5;doly=0.6;winid=-1;win=-1;nbits=3;
levels = 2^nbits;qlevels = 2.0*pi / levels; % compute quantization levels
NRounds=0;Ntrials=1;plotOrNot=1;
rot=zeros(Nsub,Np);msl_ratio=zeros(1,levels);
%% Plot Selection
show_cfg=0; %Show Configurations
%%
yyy=([0:(Nyr-1)]-(Nyr-1)/2)*doly;xxx=([0:(Nxr-1)]-(Nxr-1)/2)*dolx;
[agrdx,agrdy]=meshgrid(xxx,yyy);
%% Super Array pattern
load('june25.mat');
tx=F.txa;ty=F.tya;rot=F.rot;
tx(:,1)=tx(:,cfg_sel);ty(:,1)=ty(:,cfg_sel);rot(:,1)=rot(:,cfg_sel);
% ty(:,2)=[2.4*ones(1,4),-2.4*ones(1,4)];tx(:,2)=2*[linspace(-3,3,4),linspace(-3,3,4)];rot(:,4)=[pi*ones(1,4),zeros(1,4)];
%% Custom configs
ty(:,2)=[1.2*ones(1,4),-1.2*ones(1,4)];tx(:,2)=[linspace(-3,3,4),linspace(-3,3,4)];rot(:,2)=[pi*ones(1,4),zeros(1,4)];
%% Array of subarray pattern
for i=1:Np
    tmp=repmat(agrdx(:),[1,8])+ones(16,1)*tx(:,i).';
    txf(:,i)=tmp(:)-mean(tmp(:));
    tmp=repmat(agrdy(:),[1,8])+ones(16,1)*ty(:,i).';
    tyf(:,i)=tmp(:)-mean(tmp(:));
end
[txf2,tyf2]=fill_aperture(txf(:,1),tyf(:,1),dolx,doly);
txf2=txf2(:)-mean(txf2);tyf2=tyf2(:)-mean(tyf2);

%% Superarray Configurations
if show_cfg
    for cfg=1:Np
        figure(12)
        subplot(2,2,cfg)
        hold on;grid on;
        % for p=1:4
        %     rectangle('Position',[tx(p,cfg)-mean(txf(:,cfg))-1,-2.4,2,4.8],'EdgeColor','r','FaceColor','c')
        % end
        plot(txf(:,cfg)-mean(txf(:,cfg)),tyf(:,cfg)-mean(tyf(:,cfg)),'k.');
        title(config(cfg));xlabel('x - \lambda units');ylabel('y - \lambda units');axis([-10, 10, -10, 10]);
        %     figure(32)
        %     subplot(2,2,cfg)
        %     for i=1:Nsub
        %         xw=abs(2*cos(rot(i,cfg)))+abs(5.04*sin(rot(i,cfg)));
        %         yw=abs(2*sin(rot(i,cfg)))+abs(5.04*cos(rot(i,cfg)));
        %         image([tx(i,cfg)+1.32*sin(rot(i,cfg))-xw/2 tx(i,cfg)+1.32*sin(rot(i,cfg))+xw/2],[ty(i,cfg)-1.32*cos(rot(i,cfg))+yw/2 ty(i,cfg)-1.32*cos(rot(i,cfg))-yw/2],imcomplement(imrotate(imcomplement(img),rot(i,cfg)*180/pi))); hold on;
        %         axis([-12, 12, -12, 15]);grid on;
        %     end
        %     scatter(tx(:,cfg),ty(:,cfg),'r.')
        %     set(gca,'YDir','normal');title(config(cfg));xlabel('x - \lambda units');ylabel('y - \lambda units')
    end
end
%NOMP params
overSamplingRate_1=4;overSamplingRate_2=4;
R_s=2;R_c=3;
P_fa=0.1;
% BPDN params
solver='sdpt3';%'sdpt3'
Ns1=128;

posac=imp_samples(Ns1,0.25);% Dictionary DoAs CircPack
[posa4x,posa4y]=ndgrid(linspace(0,1-1/16,16)-0.5,linspace(0,1-1/8,8)-0.5);
posa(1,:)=posa4x(:);
posa(2,:)=posa4y(:);
% Oversampled
[posa4x,posa4y]=ndgrid(linspace(0,1,16*overSamplingRate_1)-0.5,...
    linspace(0,1,16*overSamplingRate_2)-0.5);
maskp=4*(posa4x.^2+posa4y.^2)<1;
posa_ov(1,:)=2*posa4x(maskp);
posa_ov(2,:)=2*posa4y(maskp);
posa_ov(3,:)=sqrt(1-posa_ov(1,:).^2-posa_ov(2,:).^2);
posa(3,:)=sqrt(1-posa(1,:).^2-posa(2,:).^2);
posac(3,:)=sqrt(1-posac(1,:).^2-posac(2,:).^2);
% posac=posa;% WARNING: Test
H = convhulln(posac');
Phi_el=[1,-1,1j,-1j];

figure(8)
hold on;
scatter(posa(1,:),posa(2,:),[],posa(3,:),'r+');hold on;
scatter(posac(1,:),posac(2,:),[],posac(3,:),'bo');hold on;
rectangle('Position',[-0.0625,-0.125,0.125,0.25]);
scatter(DoAs(1,:),DoAs(2,:),[],DoAs(3,:),'k*');grid on;
scatter(posa_ov(1,:),posa_ov(2,:),[],posa_ov(3,:),'md');hold on;
legend('DFT Basis','Spherical Basis','DoAs');
title('Simulation scenario, Basis & Sample DoAs');
%% JL Lemma
JL_lemma_plt=0;
if JL_lemma_plt
    % for K_DoA=1:1
    for cfg=1:Np
            Ns=size(txf,1);
        Psi_Mat = 1/sqrt(Ns)*exp(1j*2*pi*([txf(:,cfg),tyf(:,cfg)]*posa(1:2,:)));
        for i=4:4:64
            Phi_Mat0 = Phi_el(randi(4,i,Ns,N_Phi_Mat));%randn(i,Ns)+1j*randn(i,Ns);%kron(eye(8),randn(1,16));
            %         A_c = Phi_Mat0*Psi_Mat;
            for p=1:N_DoA
                Xn=Psi_Mat(:,randi(Ns,1,2*K_DoA))*[ones(K_DoA,1);-ones(K_DoA,1)];
                Yn=squeeze(Phi_Mat0(:,:,rem(p,N_Phi_Mat)+1))*Xn;
                rat(p)=norm(Yn)^2/norm(Xn)^2/i;
            end
            m_vec(i/4)=i;
            rat_vec(:,i/4)=10*log10([max(rat);min(rat)]);
            %         CEu(cfg,i)=sqrt(trace(covc));CEv(cfg,i)=sqrt(det(covc));
        end
        figure(16)
        subplot(1,2,cfg)
        plot(m_vec,rat_vec(1,:),'rx-','linewidth',K_DoA);hold on;
        plot(m_vec,rat_vec(2,:),'bx--','linewidth',K_DoA);grid on;
        %      figure(44)
        %     semilogy(snrdr,CEu(cfg,:),calcium3{cfg});hold on;grid on;
        %      figure(45)
        %     semilogy(snrdr,CEv(cfg,:),calcium3{cfg});hold on;grid on;
        axis([2,64,-9,6]);
        xlabel('M');ylabel('dB scale');title('|\Phi\Psi u|^2/(M|\Psi u|^2)');
    end
    % end
end
%% Tangent plane isometry

%% Compressive
Nalg=2;ti=zeros(Nalg,1);
fine_BPDN=1;
CEu = zeros(Np,Nsnr,Nalg);CEv=CEu;
%% Per Subarray Compressive
N_meas_sub=M_Phi/8;Ns=128;
if per_sub==1
    for ui=1:N_Phi_Mat
        if rade
            Phi_Mat_sub = Phi_el(randi(4,N_meas_sub,16,8));%kron(eye(8),randn(1,16));
        else
            Phi_Mat_sub = randn(N_meas_sub,16,8)+1i*randn(N_meas_sub,16,8);%kron(eye(8),randn(1,16));
        end
        for uj=1:8
            Phi_Mat((N_meas_sub*(uj-1)+1):N_meas_sub*uj,(16*(uj-1)+1):16*uj,ui)=Phi_Mat_sub(:,:,uj);
        end
    end
else
    if rade
        Phi_Mat = Phi_el(randi(4,M_Phi,Ns,N_Phi_Mat));%kron(eye(8),randn(1,16));
    else
        Phi_Mat = randn(M_Phi,Ns,N_Phi_Mat)+1i*randn(M_Phi,Ns,N_Phi_Mat);
    end
end
for idx=1:N_Phi_Mat
    Phi_Mat(:,:,idx)=Phi_Mat(:,:,ceil(idx/4))./kron(ones(M_Phi,1),sqrt(sum(abs(Phi_Mat(:,:,ceil(idx/4))).^2,1)));
end
%% Phi for Dense
Ns=length(txf2);
    if rade
        Phi_Mat2 = Phi_el(randi(4,M_Phi,Ns,N_Phi_Mat));%kron(eye(8),randn(1,16));
    else
        Phi_Mat2 = randn(M_Phi,Ns,N_Phi_Mat)+1i*randn(M_Phi,Ns,N_Phi_Mat);
    end
for idx=1:N_Phi_Mat
    Phi_Mat2(:,:,idx)=Phi_Mat2(:,:,ceil(idx/4))./kron(ones(M_Phi,1),sqrt(sum(abs(Phi_Mat2(:,:,ceil(idx/4))).^2,1)));
end
%% CRB
% for cfg=1:Np
%     for idx=1:N_Phi_Mat
%         [oCRBi(cfg,:,idx),oCRBu(cfg,:,idx),oCRBv(cfg,:,idx),oCRBCSi(cfg,:,idx),oCRBCSu(cfg,:,idx),...
%             oCRBCSv(cfg,:,idx)]=crb_intel(txf(:,cfg),tyf(:,cfg),snrdr,Navg,squeeze(Phi_Mat(:,:,idx)));
%     end
% end
% CRBi=sqrt(mean(oCRBi.^2,3));
% CRBu=sqrt(mean(oCRBu.^2,3));
% CRBv=sqrt(mean(oCRBv.^2,3));
% CRBCSi=sqrt(mean(oCRBCSi.^2,3));
% CRBCSu=sqrt(mean(oCRBCSu.^2,3));
% CRBCSv=sqrt(mean(oCRBCSv.^2,3));

%% Algorithms
% indx=randi(N_DoA,N_DoA2,K_DoA);
indx=(1:N_DoA2).'*ones(1,K_DoA);
% g_array=(ones(N_DoA2,1)*linspace(1/K_DoA,1,K_DoA)).*exp(2*pi*1i*rand(N_DoA2,K_DoA));
% if K_DoA==1
%     g_array=(ones(N_DoA,1)).*exp(2*pi*1i*rand(N_DoA,K_DoA));
% else
%     g_array=(ones(N_DoA,1)*[sqrt(0.1),1]).*exp(2*pi*1i*rand(N_DoA,K_DoA));
% end
% g_array=(1+rand(N_DoA2,K_DoA)).*exp(2*pi*1i*rand(N_DoA2,K_DoA));
DoAmin=-0.2; %Minimum DoA separation WARNING
mCRBu_more=zeros(Np,Nsnr,K_DoA,N_DoA,N_Phi_Mat);
mCRBCSu_more=zeros(Np,Nsnr,K_DoA,N_DoA,N_Phi_Mat);
for p=1:N_DoA2
    ua0(p,:)=DoAs(1,indx(p,:));va0(p,:)=DoAs(2,indx(p,:));
    if K_DoA>1,     ua0(p,1)=0;va0(p,1)=0;end% WARNING: Fixing 1st target to origin
    dist0(p,:)=pdist([ua0(p,:).',va0(p,:).']);
    if min(dist0(p,:))<DoAmin
        while min(dist0(p,:)) < DoAmin
            ua0(p,:)=DoAs(1,randi(N_DoA2,1,K_DoA));va0(p,:)=DoAs(2,randi(N_DoA2,1,K_DoA));
            dist0(p,:)=pdist([ua0(p,:).',va0(p,:).']);
        end
    end
end
ua=repmat(ua0,[Navg,1]);
va=repmat(va0,[Navg,1]);
for p=1:N_DoA
    for cfg=1:Np
        for phi_id=1:N_Phi_Mat
            p2=(p-1)*N_Phi_Mat+phi_id;
                [mCRBi(cfg,:,:,p2),mCRBu(cfg,:,:,p2),mCRBv(cfg,:,:,p2),mCRBCSi(cfg,:,:,p2),mCRBCSu(cfg,:,:,p2),mCRBCSv(cfg,:,:,p2)]=...
                    crb_gen_intf(txf(:,cfg),tyf(:,cfg),snrdr,N,squeeze(Phi_Mat(:,:,phi_id)),g_array(p,:),ua(p,:),va(p,:));
            
            mCRBu_more(cfg,:,:,p,phi_id)=squeeze(mCRBu(cfg,:,:,p2));
            mCRBCSu_more(cfg,:,:,p,phi_id)=squeeze(mCRBCSu(cfg,:,:,p2));
        end
    end
    
end
%% Plot CRB for individual DoA samples
if 0 % Plots across SNRs, cfgs,K_DOAs
    kd=2;% Only for target 2
    figure(46)
    for cfg=1:Np
        subplot(K_DoA,Np,cfg)
        scatter(ua(:,kd),va(:,kd),[],5*log10(mean(squeeze(abs(mCRBu_more(cfg,end,kd,:,:)).^2),2)),'filled','SizeData',100);colorbar;
        %     caxis(10*log10([min(vec(mCRBu_more(:,kd,:,:))), max(vec(mCRBu_more(:,kd,:,:)))]));
        title(['CRB - ',config{cfg}]);
        subplot(K_DoA,Np,Np+cfg)
        scatter(ua(:,kd),va(:,kd),[],5*log10(mean(squeeze(abs(mCRBCSu_more(cfg,end,kd,:,:)).^2),2)),'filled','SizeData',100);colorbar;
        %     caxis(10*log10([min(vec(mCRBCSu_more(:,kd,:,:))), max(vec(mCRBCSu_more(:,kd,:,:)))]));
        title(['CRB (CS) - ',config{cfg}]);
    end
    figure(47)
    for si=1:Nsnr
        for cfg=1:Np
            for kd=K_DoA:K_DoA
                subplot(Np*Nsnr,Nsnr,si+(cfg-1)*Nsnr)
                scatter(ua(:,kd),va(:,kd),[],5*log10(mean(squeeze(abs(mCRBu_more(cfg,si,kd,:,:)).^2),2)),'filled','SizeData',50);colorbar;
                caxis(10*log10([min(vec(mCRBu_more(:,si,kd,:,:))), max(vec(mCRBu_more(:,si,kd,:,:)))]))
                if cfg==1,     title(['Target2 SNR=',num2str(snrdr(si)),' dB']);end
                if si==1, ylabel([config{cfg},' CRB']); end
                subplot(Np*Nsnr,Nsnr,si+(cfg+Np-1)*Nsnr)
                scatter(ua(:,kd),va(:,kd),[],5*log10(mean(squeeze(abs(mCRBCSu_more(cfg,si,kd,:,:)).^2),2)),'filled','SizeData',50);colorbar;
                
                caxis(10*log10([min(vec(mCRBCSu_more(:,si,kd,:,:))), max(vec(mCRBCSu_more(:,si,kd,:,:)))]));
                if si==1, ylabel([config{cfg},' CRB(CS)']); end
            end
        end
    end
end
%% Aggregate CRB over all DoAs
CRBu2(:,:,:,1)=5*log10(mean(abs(mCRBu).^2,4));
CRBCSu2(:,:,:,1)=5*log10(mean(abs(mCRBCSu).^2,4));
CRBu2(:,:,:,2)=5*log10(mean(abs(mCRBu(:,:,:,n_far:end)).^2,4));
CRBCSu2(:,:,:,2)=5*log10(mean(abs(mCRBCSu(:,:,:,n_far:end)).^2,4));
figure(44)
for alg=1:Nalg
    subplot(1,Nalg,alg)
    for cfg=1:Np
        hold on;
        for jd=1:K_DoA
            if CS
                semilogy(snrdr,CRBCSu2(cfg,:,jd),calcium1{cfg},'Linewidth',jd);hold on;
            else
                semilogy(snrdr,CRBu2(cfg,:,jd),calcium1{cfg},'Linewidth',jd);hold on;
            end
        end
    end
    title(['RMSE vs SNR ',Alg_name{alg}]);xlabel('SNR(dB)');ylabel('Error','Interpreter','Latex');
    axis([min(snrdr),max(snrdr),-35,0]);grid on;
end
%%
figure(45)
for alg=1:Nalg
    for i_row=1:2
        subplot(2,Nalg,(i_row-1)*Nalg+alg)
        for cfg=1:Np
            hold on;
            for jd=1:K_DoA
                if CS
                    semilogy(snrdr,CRBCSu2(cfg,:,jd,i_row),calcium1{cfg},'Linewidth',jd);hold on;
                else
                    semilogy(snrdr,CRBu2(cfg,:,jd,i_row),calcium1{cfg},'Linewidth',jd);hold on;
                end
            end
        end
        axis([min(snrdr),max(snrdr),-40,0]);grid on;
    end
    title(['RMSE vs SNR ',Alg_name{alg}]);xlabel('SNR(dB)');ylabel('Error','Interpreter','Latex');
end
%%
figure(63)
for alg=1:Nalg
    subplot(1,Nalg,alg)
    for cfg=1:Np
        hold on;
        for jd=1:K_DoA
            for isnr=1:Nsnr
                if uniform_DoA
                    for ir=1:N_rings
                        crb_vs_doa(ir)= 5*log10(mean(abs(mCRBu(cfg,isnr,jd,ring_start_ind(ir):ring_start_ind(ir+1)-1)).^2));
                    end
                else
                    crb_vs_doa=5*log10(mean(reshape(abs(mCRBu(cfg,isnr,jd,1:N_DoA2)),[N_DoA2/N_rings,N_rings]).^2,1));
                end
                plot(doa_sep,crb_vs_doa,calcium5{cfg},'Linestyle','--','Linewidth',jd);
            end
        end
    end
    title(['RMSE vs DoA separation ',Alg_name{alg}]);xlabel('DoA separation (\Delta U^2+\Delta V^2)');ylabel('RMSE','Interpreter','Latex');
    axis([min(doa_sep),max(doa_sep),-35,-5]);grid on;
end
%%
fab_noise1=randn(M_Phi,N_DoA)+1i*randn(M_Phi,N_DoA);
fab_noise2=randn(Ns,N_DoA)+1i*randn(Ns,N_DoA);
cov_all=zeros(Nsnr,N_DoA2,N_Phi_Mat,Np,Nalg,K_DoA,Navg);
for cfg=1:Np
    
        Ns=size(txf,1);
        txfcfg=txf(:,cfg);tyfcfg=tyf(:,cfg);
        Phi_Matcfg=Phi_Mat;
    
    Psi_Mat = 1/sqrt(Ns)*exp(1j*2*pi*([txfcfg,tyfcfg]*posac(1:2,:)));
    if CS
        tau_scale_NOMP=log(M_Phi)-log(log(1/(1-P_fa)));
        tau_scale_Lasso=(1+1/log(M_Phi))*sqrt(log(M_Phi*4*pi*log(M_Phi)));
    else
        tau_scale_Lasso=(1+1/log(Ns))*sqrt(log(Ns*4*pi*log(Ns)));
        tau_scale_NOMP=log(Ns)-log(log(1/(1-P_fa)));
    end
    covcf=zeros(2,2,Nsnr*5,Nalg);bin_ct=zeros(Nsnr*5,1);
    covcf2=zeros(2,2,K_DoA,Nsnr,Nalg);
    covcf2var=zeros(2,2,K_DoA,Nsnr,Nalg);
    covcf2var(1,1,:,:,:)=realmax*ones(1,1,K_DoA,Nsnr,Nalg);
    fprintf('\n Running cfg=%d, SNR= ',cfg);
    for si=1:Nsnr
        covc=zeros(2,2,Nalg);
            sigma=db2mag(-snrdr(si));%/sqrt(Ns)
       
        count=1;
        h=waitbar(0,sprintf('cfg=%d, SNR= %d. Running...',cfg,snrdr(si)));
        for p=1:N_DoA
            waitbar(p/N_DoA,h);
            g0=g_array(p,:).';
            u0=ua(p,:);v0=va(p,:);
            Xn0=zeros(Ns,1);
            for kd=1:K_DoA
                    Xn0=Xn0+g0(kd)*(exp(1j*2*pi*(txfcfg*u0(kd)+tyfcfg*v0(kd))));
            end
            for phi_id=1:N_Phi_Mat
                p2=(p-1)*N_Phi_Mat+phi_id;
                Yn0=squeeze(Phi_Matcfg(:,:,phi_id))*Xn0;
                if 0
                    Yn=Yn0+sigma/sqrt(2)*sqrt(Ns/M_Phi)*fab_noise1(:,p);
                    Xn=Xn0+sigma/sqrt(2)*fab_noise2(:,p);
                else
                    Yn=Yn0+sigma/sqrt(2)*sqrt(Ns/M_Phi)*(randn(size(Yn0))+1i*randn(size(Yn0)));
                    Xn=Xn0+sigma/sqrt(2)*(randn(size(Xn0))+1i*randn(size(Xn0)));
                end
%                 if p<2, fprintf('\n SNR0=%2.2f, SNR=%2.2f, SNR_CS=%2.2f Alg:',snrdr(si),snr(Xn0,Xn-Xn0),snr(Yn0,Yn-Yn0)); end
                
                %% NOMP
                tic;
%                 fprintf('1');
                if CS
                    [omega_est, g_svd, ~] =NOMP_2dE(Yn,squeeze(Phi_Matcfg(:,:,phi_id)),txfcfg,tyfcfg,K_DoA,R_s,R_c,posa_ov);
                else
                    [omega_est, g_svd, ~] =NOMP_2dE(Xn,eye(Ns),txfcfg,tyfcfg,K_DoA,R_s,R_c,posa_ov);
                end
                u_val=omega_est(1:K_DoA,:).'/2/pi;alg=1;
                [cov_all(si,rem(p-1,N_DoA2)+1,phi_id,cfg,alg,:,ceil(p/N_DoA2)),covc(:,:,alg),covcf2(:,:,:,si,alg),covcf2var(:,:,:,si,alg)]=...
                    compute_metrics(covc(:,:,alg),covcf2(:,:,:,si,alg),squeeze(covcf2var(:,:,:,si,alg)),u_val,K_DoA,u0,v0);
%                 fprintf('2');
                o_s(cfg,si,1,count)=size(omega_est,1);
                if CS
                    [omega_est, g_svd, ~] =NOMP_2dF(Yn,squeeze(Phi_Matcfg(:,:,phi_id)),txfcfg,tyfcfg,tau_scale_NOMP*M_Phi*(sigma*sqrt(Ns/M_Phi))^2,R_s,R_c,posa_ov);
                    tau=1.1*M_Phi*(sigma*sqrt(Ns/M_Phi))^2;
                else
                    [omega_est, g_svd, ~] =NOMP_2dF(Xn,eye(Ns),txfcfg,tyfcfg,tau_scale_NOMP*Ns*(sigma)^2,R_s,R_c,posa_ov);
                    tau=1.1*Ns*(sigma)^2;
                end
                o_s(cfg,si,2,count)=size(omega_est,1);
                
                u_val2=omega_est(:,:).'/2/pi;alg=2;
                [cov_all(si,rem(p-1,N_DoA2)+1,phi_id,cfg,alg,:,ceil(p/N_DoA2)),covc(:,:,alg),covcf2(:,:,:,si,alg),covcf2var(:,:,:,si,alg)]=...
                    compute_metrics(covc(:,:,alg),covcf2(:,:,:,si,alg),squeeze(covcf2var(:,:,:,si,alg)),u_val2,K_DoA,u0,v0);
                ti(1)=ti(1)+toc;
                %% BPDN
                tic;
                if Nalg>2
%                 fprintf('3');
                % Coarse BPDN
                A_c = squeeze(Phi_Matcfg(:,:,phi_id))*Psi_Mat;
                if CS
                    [u_hat,~]=CVXL1(A_c,Yn,sqrt(M_Phi)*tau_scale_Lasso*sigma*sqrt(Ns/M_Phi),solver);
                else
                    [u_hat,~]=CVXL1(Psi_Mat,Xn,sqrt(Ns)*tau_scale_Lasso*sigma,solver);
                end
                active=abs(u_hat)>0.05*max(abs(u_hat));%was 0.05
                posa3=posac;Psi_Mat2=Psi_Mat;H3=H;
                active3=active;u_hat3=u_hat;
                % Fine BPDN
                for ref=1:Nref
                    %                         H = convhulln(posa3');
                    %Plotting
                    %                         figure(66)
                    %                         subplot(1,Nref,ref)
                    %%
                    N_H3=size(H3,1);H4=H3;del_indx=[];
                    for i=1:N_H3
                        if any(active3(H3(i,:)))
                            v = sum(posa3(:,H3(i,:)),2);
                            v = v / norm(v,2);
                            posa3(:,end+1) = v;
                            Psi_Mat2(:,end+1)= 1/sqrt(Ns)*exp(1j*2*pi*(txfcfg*v(1)+tyfcfg*v(2)));
                            posa_N=size(posa3,2);
                            %                                 current_N_H4=size(H4,1);
                            %                                 H4(current_N_H4+1,:)=[H3(i,1:2),posa_N];% Add 3 facets to hull
                            %                                 H4(current_N_H4+2,:)=[H3(i,2:3),posa_N];
                            %                                 H4(current_N_H4+3,:)=[H3(i,1:2:3),posa_N];
                            %                                 del_indx=[del_indx,i];
                        end
                        %                             patch(posa3(1,H3(i,:)),posa3(2,H3(i,:)),posa3(3,H3(i,:)),mag2db(abs(u_hat(H3(i,:)))),'FaceColor','interp');hold on;
                    end
                    %                         H4(del_indx,:)=[];%Delete current facet of convex hull
                    %                         H3=H4;
                    A_cf = squeeze(Phi_Matcfg(:,:,phi_id))*Psi_Mat2;
                    if CS
                        u_hat3=CVXL1(A_cf,Yn,sqrt(M_Phi)*tau_scale_Lasso*sigma*sqrt(Ns/M_Phi),solver);
                    else
                        u_hat3=CVXL1(Psi_Mat2,Xn,sqrt(Ns)*tau_scale_Lasso*sigma,solver);
                    end
                    active3=abs(u_hat3)>0.05*max(abs(u_hat3));%was 0.05
                    %                         ind_select=[1:Ns1,Ns1+find(active(Ns1+1:end)).'];
                    %                         posa3=posa3(:,ind_select);Psi_Mat2=Psi_Mat2(:,ind_select);
                    %                         active=active(ind_select);
                    if isnan(u_hat), break; else active=active3;u_hat=u_hat3; end
                    if nnz(active)>256, break; end
                    H3=convhulln(posa3');
                end
                o_s(cfg,si,3,count)=nnz(active);
                
                %                     return;
                %                     active=abs(u_hat)>0.05*max(abs(u_hat));
                omegaList=posa3(1:2,active).'*2*pi;
                %                     gainList=u_hat(ind_select(active));
                gainList=u_hat(active);
                if size(gainList,1)>128
                    [~,id]=sort(abs(gainList),'descend');
                    gainList=gainList(id(1:100));
                    omegaList=omegaList(id(1:100),1:2);
                end
                
                %%
                u_val3=omegaList.'/2/pi;alg=3;
                [cov_all(si,rem(p-1,N_DoA2)+1,phi_id,cfg,alg,:,ceil(p/N_DoA2)),covc(:,:,alg),covcf2(:,:,:,si,alg),covcf2var(:,:,:,si,alg)]=...
                    compute_metrics(covc(:,:,alg),covcf2(:,:,:,si,alg),squeeze(covcf2var(:,:,:,si,alg)),u_val3,K_DoA,u0,v0);
                end
                                %% NLasso
                if Nalg>3
%                 fprintf('4');
                %% Clustering
                
                if size(omegaList,1)>2
                    K_BIC=2;
                    [~,omegaListc]=kmeans(omegaList,K_BIC);
                    gainListc=abs(omegaListc(:,1));% Just a placeholder
                else
                    omegaListc=omegaList;gainListc=gainList;
                end
                if CS
                    [omega_NLasso, g_NLasso, residue_NLasso]=NLasso(Yn,Phi_Matcfg(:,:,phi_id),txfcfg,tyfcfg,tau_scale_Lasso,R_c,omegaListc,gainListc);
                else
                    [omega_NLasso, g_NLasso, residue_NLasso]=NLasso(Xn,eye(Ns),txfcfg,tyfcfg,tau_scale_Lasso,R_c,omegaListc,gainListc);
                end
                u_val4=omega_NLasso.'/2/pi;alg=4;
                [cov_all(si,rem(p-1,N_DoA2)+1,phi_id,cfg,alg,:,ceil(p/N_DoA2)),covc(:,:,alg),covcf2(:,:,:,si,alg),covcf2var(:,:,:,si,alg)]=...
                    compute_metrics(covc(:,:,alg),covcf2(:,:,:,si,alg),squeeze(covcf2var(:,:,:,si,alg)),u_val4,K_DoA,u0,v0);
                %% Choosing peak on either side of decision boundary
                indxtarget=([u0(2),v0(2)]*posa3(1:2,:)-(u0(2)^2+v0(2)^2)/2)>0;
                [~,u_s]=max(abs(u_hat(~indxtarget)));
                posns=posa3(1:2,~indxtarget);
                u_val(:,1)=posns(1:2,u_s);
                [~,u_t]=max(abs(u_hat(indxtarget)));
                posnt=posa3(1:2,indxtarget);
                u_val(:,2)=posnt(1:2,u_t);
                o_s(cfg,si,4,count)=size(omega_NLasso,1);
                %% Choosing largest 2 values
                %                     [~,u_s]=max(abs(u_hat));
                %                     u_val(:,1)=posa3(1:2,u_s);
                %                     u_hat(u_s)=-realmax;
                %                     [~,u_t]=max(abs(u_hat));
                %                     u_val(:,2)=posa3(1:2,u_t);
                %                 hold on;plot3(u_val(1,:),u_val(2,:),[50,50],'rs');
                end
                %%
                ti(2)=ti(2)+toc;
                count=count+1;
            end
        end
        delete(h);
        covc=covc/N_DoA/N_Phi_Mat;
        for alg=1:Nalg
            CEu(cfg,si,alg)=5*log10(trace(covc(:,:,alg)));CEv(cfg,si,alg)=5*log10(det(covc(:,:,alg)));
        end
    end
    figure(44)
    for alg=1:Nalg
        subplot(1,Nalg,alg)
        plot(snrdr,CEu(cfg,:,alg),calcium3{cfg});hold on;grid on;
    end

    [cov_avg(1:Nsnr,cfg,1:Nalg,1:K_DoA),cov_min(1:Nsnr,cfg,1:Nalg,1:K_DoA),cov_max(1:Nsnr,cfg,1:Nalg,1:K_DoA)]=...
        collect_stats(cov_all,1:size(cov_all,2),cfg,Nsnr,Nalg,K_DoA);

    [cov_avg2(1:Nsnr,cfg,1:Nalg,1:K_DoA),cov_min2(1:Nsnr,cfg,1:Nalg,1:K_DoA),cov_max2(1:Nsnr,cfg,1:Nalg,1:K_DoA)]=...
        collect_stats(cov_all,n_far:size(cov_all,2),cfg,Nsnr,Nalg,K_DoA);
    % [cov_avg3,cov_min3,cov_max3]=collect_stats(cov_all,n_far2:size(cov_all,2));
    cov_avg_N=zeros(Nsnr,cfg,Nalg,K_DoA,N_rings);cov_max_N=cov_avg_N;
    for ir=1:N_rings% Collect RMSE vs DoA separation
        if uniform_DoA
                    [cov_avg_N(1:Nsnr,cfg,1:Nalg,1:K_DoA,ir),~,cov_max_N(1:Nsnr,cfg,1:Nalg,1:K_DoA,ir)]=...
            collect_stats(cov_all,ring_start_ind(ir):ring_start_ind(ir+1)-1,cfg,Nsnr,Nalg,K_DoA);
        else
        [cov_avg_N(1:Nsnr,cfg,1:Nalg,1:K_DoA,ir),~,cov_max_N(1:Nsnr,cfg,1:Nalg,1:K_DoA,ir)]=...
            collect_stats(cov_all,N_DoA2/N_rings*(ir-1)+1:N_DoA2/N_rings*(ir),cfg,Nsnr,Nalg,K_DoA);
        end
    end
    for jd=1:K_DoA
        for alg=1:Nalg
            figure(44)
            subplot(1,Nalg,alg)
            errorbar(snrdr,5*log10(squeeze(covcf2(1,1,jd,:,alg)+covcf2(2,2,jd,:,alg))/N_DoA/N_Phi_Mat),...
                5*log10(squeeze(covcf2(1,1,jd,:,alg)+covcf2(2,2,jd,:,alg))/N_DoA/N_Phi_Mat)-5*log10(squeeze(covcf2var(1,1,jd,:,alg))),...
                5*log10(squeeze(covcf2var(2,2,jd,:,alg)))-5*log10(squeeze(covcf2(1,1,jd,:,alg)+covcf2(2,2,jd,:,alg))/N_DoA/N_Phi_Mat),...
                calcium5{cfg},'Linewidth',jd);
            figure(45)
            subplot(2,Nalg,alg)
            errorbar(snrdr,cov_avg(:,cfg,alg,jd),zeros(size(snrdr)),cov_max(:,cfg,alg,jd),calcium5{cfg},'Linewidth',jd);
            subplot(2,Nalg,Nalg+alg)
            errorbar(snrdr,cov_avg2(:,cfg,alg,jd),zeros(size(snrdr)),cov_max2(:,cfg,alg,jd),calcium5{cfg},'Linewidth',jd);
            figure(63)
            subplot(1,Nalg,alg)
            for isnr=1:Nsnr
            plot(doa_sep,squeeze(cov_avg_N(isnr,cfg,alg,jd,:)),calcium5{cfg},'Linewidth',jd);hold on;
            end
        end
    end
end
text(-8,-45,[dir_name,batch_name],'Interpreter','none');
%%
figure(44)
subplot(1,Nalg,1)
legend(vec([strcat(config,' -CRB(CS)1');strcat(config,' -CRB(CS)2');strcat(config,' -RMSEavg');strcat(config,' -RMSE1');strcat(config,' -RMSE2')]),'location','northeast');hold on;
subplot(1,Nalg,2)
legend(vec([strcat(config,' -CRB1');strcat(config,' -CRB2');strcat(config,' -RMSEavg');strcat(config,' -RMSE1');strcat(config,' -RMSE2')]),'location','northeast');hold on;
%% Plot CRB for individual DoA samples
%%
errorplots=0;
if errorplots
    cov_all2(1:Nsnr,1:N_DoA2,1:Np,1:Nalg)=5*log10(squeeze(mean(mean(cov_all(:,:,:,:,:,2,:),3),7)));
    normz1=max(vec(cov_all2(end,:,:,1)));
    normz2=max(vec(cov_all2(end,:,:,2)));
    rnge=(1:N_DoA2);
    figure(48)
    for cfg=1:Np
        for kd=K_DoA:K_DoA
            subplot(Np,2,cfg)
            scatter(ua(rnge,kd),va(rnge,kd),[],squeeze(cov_all2(end,:,cfg,2)).','filled');colorbar;
            caxis([min(vec(cov_all2(end,rnge,:,2))), max(vec(cov_all2(end,rnge,:,2)))])
            title(['Target2 RMSE SNR= ',num2str(snrdr(end)),config{cfg},'(M=N)']);
            
            subplot(Np,2,Np+cfg)
            scatter(ua(rnge,kd),va(rnge,kd),[],squeeze(cov_all2(end,:,cfg,1)).','filled');colorbar;
            caxis([min(vec(cov_all2(end,rnge,:,1))), max(vec(cov_all2(end,rnge,:,1)))])
            title(['Target2 RMSE SNR= ',num2str(snrdr(end)),config{cfg},'(M=',num2str(M_Phi)]);
        end
    end
    figure(49)
    for si=1:Nsnr
        for cfg=1:Np
            for kd=K_DoA:K_DoA
                subplot(4,Nsnr,si+(cfg-1)*Nsnr)
                scatter(ua(rnge,kd),va(rnge,kd),[],squeeze(cov_all2(si,:,cfg,2)).','filled','SizeData',60);colorbar;
                caxis([min(vec(cov_all2(si,:,:,2))), max(vec(cov_all2(si,:,:,2)))])
                if cfg==1,     title(['Target2 SNR=',num2str(snrdr(si)),' dB']);end
                if si==1, ylabel([config{cfg},' (M=N)']); end
                subplot(4,Nsnr,si+(cfg+Np-1)*Nsnr)
                scatter(ua(rnge,kd),va(rnge,kd),[],squeeze(cov_all2(si,:,cfg,1)).','filled','SizeData',60);colorbar;
                caxis([min(vec(cov_all2(si,:,:,1))), max(vec(cov_all2(si,:,:,1)))])
                if si==1, ylabel([config{cfg},' (M=',num2str(M_Phi),')']); end
            end
        end
    end
    %%
    figure(50)
    for si=1:N_Phi_Mat
        for cfg=1:Np
            for kd=K_DoA:K_DoA
                subplot(4,N_Phi_Mat,si+(cfg-1)*N_Phi_Mat)
                scatter(ua(rnge,kd),va(rnge,kd),[],5*log10(squeeze(mean(cov_all(end,:,si,cfg,2,2,:),7))).','filled','SizeData',200);colorbar;
                caxis(5*log10([min(vec(cov_all(end,:,si,:,2,2,:))), max(vec(cov_all(end,:,si,:,2,2,:)))]));
                %            title(['RMSE (M=N)',config{cfg}]);
                
                subplot(4,N_Phi_Mat,si+(cfg+Np-1)*N_Phi_Mat)
                scatter(ua(rnge,kd),va(rnge,kd),[],5*log10(squeeze(mean(cov_all(end,:,si,cfg,1,2,:),7))).','filled','SizeData',200);colorbar;
                caxis(5*log10([min(vec(cov_all(end,:,si,:,1,2,:))), max(vec(cov_all(end,:,si,:,1,2,:)))]));
                %                    title(['RMSE (M=',num2str(M_Phi),')',config{cfg}]);
            end
        end
    end
    %%
    cov_all1(1:Nsnr,1:N_DoA2,1:Np,1:Nalg)=5*log10(squeeze(mean(mean(cov_all(:,:,:,:,:,1,:),3),7)));
    figure(51)
    for cfg=1:Np
        for kd=K_DoA:K_DoA
            subplot(Np,2,cfg)
            scatter(ua(rnge,kd),va(rnge,kd),[],squeeze(cov_all1(end,:,cfg,2)).','filled');colorbar;
            caxis([min(vec(cov_all1(end,:,:,2))), max(vec(cov_all1(end,:,:,2)))])
            title(['Target1 RMSE SNR= ',num2str(snrdr(end)),config{cfg},'(M=N)']);
            subplot(Np,2,Np+cfg)
            scatter(ua(rnge,kd),va(rnge,kd),[],squeeze(cov_all1(end,:,cfg,1)).','filled');colorbar;
            caxis([min(vec(cov_all1(end,:,:,1))), max(vec(cov_all1(end,:,:,1)))])
            title(['Target1 RMSE SNR= ',num2str(snrdr(end)),config{cfg},'(M=',num2str(M_Phi)]);
        end
    end
    %%
    figure(52)
    for si=1:Nsnr
        for cfg=1:Np
            for kd=K_DoA:K_DoA
                subplot(4,Nsnr,si+(cfg-1)*Nsnr)
                scatter(ua(rnge,kd),va(rnge,kd),[],squeeze(cov_all1(si,:,cfg,2)).','filled','SizeData',60);colorbar;
                caxis([min(vec(cov_all1(si,:,:,2))), max(vec(cov_all1(si,:,:,2)))])
                if cfg==1,     title(['Target1 SNR=',num2str(snrdr(si)),' dB']);end
                if si==1, ylabel([config{cfg},' (M=N)']); end
                subplot(4,Nsnr,si+(cfg+Np-1)*Nsnr)
                scatter(ua(rnge,kd),va(rnge,kd),[],squeeze(cov_all1(si,:,cfg,1)).','filled','SizeData',60);colorbar;
                caxis([min(vec(cov_all1(si,:,:,1))), max(vec(cov_all1(si,:,:,1)))])
                if si==1, ylabel([config{cfg},' (M=',num2str(M_Phi),')']); end
            end
        end
    end
    %%
    colr=['r','b','g','m','c','k'];
    figure(53)
    kd=K_DoA;Rc_u=0.25;
    a=2*pi*(1-sqrt(1-Rc_u))/N_DoA2;d=sqrt(a);
    M_lat=round(asin(sqrt(Rc_u))/d);
    d_lat=asin(sqrt(Rc_u))/M_lat;
    d_lon=a/d_lat;
    
    for si=1:Nsnr
        for cfg=1:Np
            idxt=1;
            for p=1:M_lat
                %             lat_a=([3*p-2,3*p-1,3*p]-0.5)*asin(sqrt(Rc_u))/M_lat;
                lat_a=(p-0.5)*asin(sqrt(Rc_u))/M_lat;
                M_lon=sum(round(2*pi*sin(lat_a)/d_lon));
                if (idxt+M_lon)>(N_DoA2+1), M_lon=N_DoA2-idxt; end
                subplot(4,Nsnr,si+(cfg-1)*Nsnr)
                [Nh,eh]=histcounts(vec(cov_all(si,idxt:idxt+M_lon-1,:,cfg,2,2,:)),'Normalization','probability');%probability
                plot((eh(1:end-1)+eh(2:end))/2,Nh,'color',colr(p));hold on;
                
                %            title(['RMSE (M=N)',config{cfg}]);50,'DisplayStyle','stairs'
                set(gca, 'XScale', 'log');set(gca, 'YScale', 'log');
                subplot(4,Nsnr,si+(cfg+Np-1)*Nsnr)
                [Nh,eh]=histcounts(vec(cov_all(si,idxt:idxt+M_lon-1,:,cfg,1,2,:)),'Normalization','probability');
                plot((eh(1:end-1)+eh(2:end))/2,Nh,'color',colr(p));hold on;
                %                    title(['RMSE (M=',num2str(M_Phi),')',config{cfg}]);
                set(gca, 'XScale', 'log');set(gca, 'YScale', 'log');
                idxta(idxt:idxt+M_lon-1)=p;
                idxt=idxt+M_lon;
                
            end
        end
    end
    %%
    % figure(61)
    % scatter(ua(1:N_DoA2,2),va(1:N_DoA2,2),[],idxta,'filled','SizeData',100);colormap('lines');
    idxt=1;
    figure(61)
    for p=1:M_lat
        lat_a=(p-0.5)*asin(sqrt(Rc_u))/M_lat;
        M_lon=sum(round(2*pi*sin(lat_a)/d_lon));
        if (idxt+M_lon)>(N_DoA2+1), M_lon=N_DoA2-idxt; end
        scatter(ua(idxt:idxt+M_lon-1,2),va(idxt:idxt+M_lon-1,2),[],colr(p),'filled','SizeData',100);
        hold on;
        idxt=idxt+M_lon;
    end
end
%%
mkr={'k--','rs-','b.-','m*-'};

figure(62)
for cfg=1:Np
    for si=1:Nsnr
        subplot(Np,Nsnr,(cfg-1)*Nsnr+si)
        for alg=1:Nalg
            [h,x]=ecdf(squeeze(o_s(cfg,si,alg,:)));
            stairs(x,h,mkr{alg});hold on;
            set(gca,'Xscale','linear');
        end
        grid on;title(['K_{DoA} CDF: SNR=',num2str(snrdr(si)),' dB, ',config{cfg}]);xlabel('K_{DoA}');ylabel('p(k<K_{DoA})');
        axis([0,10,-0.1,1.1]);
        legend(Alg_name,'location','southeast')
    end
end
text(-8,-0.2,[dir_name,batch_name],'Interpreter','none');
%% CDF's
figure(64)
for alg=1:Nalg
    
    for cfg=1:Np
        subplot(Nalg,3,3*(alg-1)+1)
        plot(snrdr,mean(squeeze(o_s(cfg,:,alg,:))==2,2),calcium0{cfg});hold on;grid on;
        title(['Correct K_{DoA}: ',Alg_name{alg}]);
        xlabel('SNR (dB)');ylabel('p(K=2)');axis([min(snrdr),max(snrdr),-0.2,1.2]);
        subplot(Nalg,3,3*(alg-1)+2)
        plot(snrdr,mean(squeeze(o_s(cfg,:,alg,:))<2,2),calcium0{cfg});hold on;
                title(['Correct K_{DoA}: ',Alg_name{alg}]);grid on;
        xlabel('SNR (dB)');ylabel('p(K<2)');axis([min(snrdr),max(snrdr),-0.2,1.2]);
        subplot(Nalg,3,3*(alg-1)+3)
        plot(snrdr,mean(squeeze(o_s(cfg,:,alg,:))>2,2),calcium0{cfg});hold on;grid on;
                title(['Correct K_{DoA}: ',Alg_name{alg}]);
        xlabel('SNR (dB)');ylabel('p(K>2)');axis([min(snrdr),max(snrdr),-0.2,1.2]);
    end
end
figure(65)
for alg=1:Nalg
    for cfg=1:Np
        subplot(Nalg,3,3*(alg-1)+1)
        for si=1:Nsnr
            p_detect=zeros(N_rings,1);
            for ir=1:N_rings% Collect RMSE vs DoA separation
                if uniform_DoA
                    for navg=1:Navg
                        p_detect(ir)=p_detect(ir)+ mean(squeeze(o_s(cfg,si,alg,(ring_start_ind(ir):ring_start_ind(ir+1)-1)+N_DoA2*(navg-1)))==2,1)/Navg;
                    end
                else
                    for navg=1:Navg
                        p_detect(ir)=p_detect(ir)+ mean(squeeze(o_s(cfg,si,alg,(N_DoA2/N_rings*(ir-1)+1:N_DoA2/N_rings*(ir))+N_DoA2*(navg-1)))==2,1);
                    end
                end
            end
                        plot(doa_sep,p_detect,calcium0{cfg},'Linewidth',si);hold on;
        end
            title(['Correct K_{DoA}: ',Alg_name{alg}]);
            xlabel('DoA separation');ylabel('p(K=2)');
        axis([0,0.4,-0.2,1.2]);grid on;
        subplot(Nalg,3,3*(alg-1)+2)
        for si=1:Nsnr
            p_detect=zeros(N_rings,1);
            for ir=1:N_rings% Collect RMSE vs DoA separation
                if uniform_DoA
                    for navg=1:Navg
                        p_detect(ir)=p_detect(ir)+ mean(squeeze(o_s(cfg,si,alg,(ring_start_ind(ir):ring_start_ind(ir+1)-1)+N_DoA2*(navg-1)))<2,1)/Navg;
                    end
                else
                    for navg=1:Navg
                        p_detect(ir)=p_detect(ir)+ mean(squeeze(o_s(cfg,si,alg,(N_DoA2/N_rings*(ir-1)+1:N_DoA2/N_rings*(ir))+N_DoA2*(navg-1)))<2,1);
                    end
                end
            end
            plot(doa_sep,p_detect,calcium0{cfg},'Linewidth',si);hold on;
        end
                    title(['Lower K_{DoA}: ',Alg_name{alg}]);
            xlabel('DoA separation');ylabel('p(K<2)');
        axis([0,0.4,-0.2,1.2]);grid on;
        subplot(Nalg,3,3*(alg-1)+3)
        for si=1:Nsnr
            p_detect=zeros(N_rings,1);
            for ir=1:N_rings% Collect RMSE vs DoA separation
                if uniform_DoA
                    for navg=1:Navg
                        p_detect(ir)=p_detect(ir)+ mean(squeeze(o_s(cfg,si,alg,(ring_start_ind(ir):ring_start_ind(ir+1)-1)+N_DoA2*(navg-1)))>2,1)/Navg;
                    end
                else
                    for navg=1:Navg
                        p_detect(ir)=p_detect(ir)+ mean(squeeze(o_s(cfg,si,alg,(N_DoA2/N_rings*(ir-1)+1:N_DoA2/N_rings*(ir))+N_DoA2*(navg-1)))>2,1);
                    end
                end
            end
                        plot(doa_sep,p_detect,calcium0{cfg},'Linewidth',si);hold on;
        end
                    title(['Higher K_{DoA}: ',Alg_name{alg}]);
            xlabel('DoA separation');ylabel('p(K>2)');
        axis([0,0.4,-0.2,1.2]);grid on;
    end
end
%%
figure(60)
bar(ti);
%%

mkdir([dir_name,batch_name]);
% saveas(figure(44),[dir_name,batch_name,'/MSEvsCRB']);
saveas(figure(8),[dir_name,batch_name,'/DoAs']);
saveas(figure(45),[dir_name,batch_name,'/MSEvsCRB_both']);
saveas(figure(62),[dir_name,batch_name,'/CDF_snrs_configs']);
saveas(figure(63),[dir_name,batch_name,'/MSEvsDoA']);
saveas(figure(64),[dir_name,batch_name,'/prob_snrs']);
saveas(figure(65),[dir_name,batch_name,'/prob_doa']);
if errorplots
    saveas(figure(46),[dir_name,batch_name,'/CRB']);
    saveas(figure(48),[dir_name,batch_name,'/MSE_last_target']);
    saveas(figure(49),[dir_name,batch_name,'/MSE_SNR_target']);
    saveas(figure(51),[dir_name,batch_name,'/MSE_last_source']);
    saveas(figure(52),[dir_name,batch_name,'/MSE_SNR_source']);
    saveas(figure(53),[dir_name,batch_name,'/MSE_hist']);
    saveas(figure(60),[dir_name,batch_name,'/runtime']);
    saveas(figure(61),[dir_name,batch_name,'/hist_plot_regions']);
end
%% CDF of pairwise DoA distances in each trial
% figure;cdfplot(dist0(:))
