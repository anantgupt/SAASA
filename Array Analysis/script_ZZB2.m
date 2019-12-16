% Use this program to plot radiation pattern for array(s) of subarray
% For configs with 8 subarray, optimize positioning
% Anant Gupta, UCSB
clc;
clear all;
close all;
% To put all the figures in tabs
set(0,'DefaultFigureWindowStyle','docked');
%%
eps =0.001;rho=1.0;
Nres=256;midi=Nres/2+1;
bws=3;%Given beamwidth(doublesided) region in degree
img = imread('rfem.png');             %# Load rfem image
norm_compare=@(in)(in-0.5*(max(in)+min(in)))/(max(in)-min(in));
configs={'Tri','crc','lin','B'};Np=3;Nsub=8;
% config={'Triangle','Circle','Linear','Benchmark'};
config={'Sparse','Compact','Dense (SNR-adjusted)'};
colr=['r','b','g','k'];
%% Subarray beam pattern
Nxr=4;Nyr=4;dolx=0.5;doly=0.6;winid=-1;win=-1;nbits=3;
levels = 2^nbits;qlevels = 2.0*pi / levels; % compute quantization levels
NRounds=0;Ntrials=1;plotOrNot=1;
rot=zeros(Nsub,Np);msl_ratio=zeros(1,levels);
U = [-Nres/2:(Nres/2)-1] ./(Nres/2);V = [-Nres/2:(Nres/2)-1] ./(Nres/2);
[SU,SV] = meshgrid(V,U);
mask2=(SU.^2+SV.^2) < 1-2/Nres;
for u=1:size(SU,1)
    for v=1:size(SU,2)
        pb(v,u)=patchbeam(SU(v,u),SV(v,u));
    end
end
pbi=20*log10(pb);
indx=find((SU.^2+SV.^2)>1);
yyy=([0:(Nyr-1)]-(Nyr-1)/2)*doly;xxx=([0:(Nxr-1)]-(Nxr-1)/2)*dolx;
[agrdx,agrdy]=meshgrid(xxx,yyy);
%% Super Array pattern
if 0 % BW,ecc
load('Sep5.mat');%apr12
tx=P.txa;ty=P.tya;rot=P.rot;
else % BW,MSLL
load('Sep18.mat');%apr12%
tx=F.txa;ty=F.tya;rot=F.rot;
end
%% Custom configs
cfg_sel=2;
date_str=['jun25_',num2str(cfg_sel)];
tx(:,1)=tx(:,cfg_sel);ty(:,1)=ty(:,cfg_sel);rot(:,1)=rot(:,cfg_sel);
ty(:,2)=[1.2*ones(1,4),-1.2*ones(1,4)];tx(:,2)=[linspace(-3,3,4),linspace(-3,3,4)];rot(:,2)=[pi*ones(1,4),zeros(1,4)];
for i=1:Np
    tmp=repmat(tx(:,i),[1,16])+ones(8,1)*agrdx(:).';
    txf(:,i)=tmp(:)-mean(tmp(:));% Center the array
    tmp=repmat(ty(:,i),[1,16])+ones(8,1)*agrdy(:).';
    tyf(:,i)=tmp(:)-mean(tmp(:));% Center the array 
end
% txf2=txf;tyf2=tyf;
[txf2,tyf2]=fill_aperture(txf(:,1),tyf(:,1),dolx,doly);
txf2=txf2-mean(txf2);tyf2=tyf2-mean(tyf2);%Centering
%% Superarray Configurations
for cfg=1:Np
    figure(12)
    subplot(1,3,cfg)
    hold on;grid on;
for p=1:4
    rectangle('Position',[tx(p,cfg)-mean(txf(:,cfg))-1,-2.4,2,4.8],'EdgeColor','r','FaceColor','c')
end
    plot(txf(:,cfg)-mean(txf(:,cfg)),tyf(:,cfg)-mean(tyf(:,cfg)),'k.');
    title(config(cfg));xlabel('x - \lambda units');ylabel('y - \lambda units');axis([-10, 10, -10, 10]);
    figure(32)
    subplot(1,3,cfg)
    for i=1:Nsub
        xw=abs(2*cos(rot(i,cfg)))+abs(5.04*sin(rot(i,cfg)));
        yw=abs(2*sin(rot(i,cfg)))+abs(5.04*cos(rot(i,cfg)));
        image([tx(i,cfg)+1.32*sin(rot(i,cfg))-xw/2 tx(i,cfg)+1.32*sin(rot(i,cfg))+xw/2],[ty(i,cfg)-1.32*cos(rot(i,cfg))+yw/2 ty(i,cfg)-1.32*cos(rot(i,cfg))-yw/2],imcomplement(imrotate(imcomplement(img),rot(i,cfg)*180/pi))); hold on;
        axis([-12, 12, -12, 15]);grid on;
    end
    scatter(tx(:,cfg),ty(:,cfg),'r.')
    set(gca,'YDir','normal');title(config(cfg));xlabel('x - \lambda units');ylabel('y - \lambda units')
end
%% ZZB
calcium={'r--v','b--','g--+','k--'};calcium2={'r-v','b-','g-+','k-'};calcium3={'r-.v','b-.','g-.+','k-.'};calcium4={'r:x','b:x','g:+','k-'};
sb=-25;se=10;
snrdr=[-14:1:-9,-8.6:0.2:-6,-5.5:0.5:-4,0];%[sb:4:-13,-12:1:-6,-4:4:0];
% snrdr=[-20:10:0];
N=1;%was 20
Np=3;Ns=64;pb_case=0;
%     ha=2*(0:0.05:1).^2;
%     vva=[-(-1:0.2:0),(0:0.2:1).^2];
load('Alut2.mat');% Contains preset ha,vva,A

for cfg=1:Np
    figure(40)
    subplot(1,3,cfg)
    if cfg<3
    [CRBi(cfg,:),CRBu(cfg,:),CRBv(cfg,:),ZZBi(cfg,:),ZZBu(cfg,:),ZZBv(cfg,:),c_CRB]=crb_intel(txf(:,cfg),tyf(:,cfg),snrdr,N);
    else
        [CRBi(cfg,:),CRBu(cfg,:),CRBv(cfg,:),ZZBi(cfg,:),ZZBu(cfg,:),ZZBv(cfg,:),c_CRB]=crb_intel(txf2,tyf2,snrdr-10*log10(length(txf2)/length(txf(:,1))),N);
end
end
%%
    for cfg=1:Np
       figure(43)
    subplot(2,1,1)
    plot(snrdr,10*log10(CRBu(cfg,:)),calcium{cfg});hold on;title('RMSE Max');xlabel('SNR(dB)');ylabel('RMSE','Interpreter','Latex');%axis([sb,se,1e-4,1]);
    subplot(2,1,2)
    plot(snrdr,10*log10(CRBv(cfg,:)),calcium{cfg});hold on;title('RMSE Min');xlabel('SNR(dB)');ylabel('RMSE','Interpreter','Latex');%axis([sb,se,1e-4,1]);
    figure(44)
    plot(snrdr,10*log10(CRBi(cfg,:)),calcium{cfg});hold on;title('Error Area');xlabel('SNR(dB)');ylabel('Error','Interpreter','Latex');
figure(45)
%     semilogy(snrdr,sqrt(CRBu(cfg,:).^2+CRBv(cfg,:).^2),calcium{cfg});hold on;title('Error Area');xlabel('SNR(dB)');ylabel('Error','Interpreter','Latex');
    plot(snrdr,5*log10(CRBu(cfg,:).^2+CRBv(cfg,:).^2),calcium{cfg});hold on;title('Overall Error');xlabel('SNR(dB)');ylabel('Error','Interpreter','Latex');

    end
%%
% figure(42)
%     subplot(1,3,3)
% legend(configs,'location','southwest')
% subplot(1,3,1)
% legend(configs,'location','southwest')
% subplot(1,3,2)
% legend(configs,'location','southwest')
% figure(43)
% subplot(1,2,1)
% legend(configs,'location','southwest')
% subplot(1,2,2)
% legend(configs,'location','southwest')
%% BB from 
% fprintf('\n BB:');
% %% Find common Test points
% u0=0;v0=0;peaksid=[];
% for cfg=1:2
%     pattern=mag2db(abs( bpmra(txf(:,cfg),tyf(:,cfg),SU,SV,u0,v0,-1,indx) ));
% mask=true(3);mask(5)=0;Ng=size(pattern,2);midi=Ng/2+1;
% mask2=(SU.^2+SV.^2) < 0.25;%1-2/N;% Removes unnecessary peaks found on boundary
% p=pattern.*(pattern>ordfilt2(pattern,8,mask)).*mask2;%.*(((U-umax).^2+(V-vmax).^2)>amin);
%     p(isnan(p))=0;% NaN were encountered for + config bcoz -Inf.
%     [peaks0,id2t]=sort(p(:),'descend');
%     peakid0=id2t(1:nnz(peaks0));
%     peaksid=union(peakid0,peaksid);
% %     figure(34)
% %     subplot(1,2,cfg)
% %     surf(SU,SV,max(pattern(:,:),0),'linestyle','none');view(2);%colorbar;%colormap(flipud(colormap));
% % hold on;scatter3(SU(peakid0),SV(peakid0),max(max(pattern))*ones(size(peakid0)),'rx');
% %             axis([-0.5, 0.5, -0.5, 0.5]);
% end
% if 0 % random test points
%     N_bb=51;
%     randindx=find(mask2>0);selindx=randi(length(randindx),N_bb-1,1);size(randindx(selindx))
%     peaksid=randindx(selindx);
% end
% if 4 % UNiformly chosen testpoints
% posa=imp_samples(Ns,0.25);
% UVd(:,1)=posa(1,:).';
% UVd(:,2)=posa(2,:).';
% else
%     UVd=[SU(peaksid),SV(peaksid)];
% end
% for cfg=1:Np
%             pattern=mag2db(abs( bpmra(txf(:,cfg),tyf(:,cfg),SU,SV,u0,v0,-1,indx) ));
% figure(34)
%     subplot(1,2,cfg)
%     surf(SU,SV,max(pattern(:,:),0),'linestyle','none');view(2);%colorbar;%colormap(flipud(colormap));
% hold on;scatter3(UVd(:,1),UVd(:,2),max(max(pattern))*ones(size(UVd(:,1))),'rx');
%             axis([-0.5, 0.5, -0.5, 0.5]);
%         end
%     
%     for i=1:1
%         figure(40)
%         for cfg=1:Np
%             [BB1(cfg,:),BBu1(cfg,:),BBv1(cfg,:)]=bb_intel(txf(:,cfg),tyf(:,cfg),snrdr,N,UVd);fprintf(' %d',cfg);
%         end
%         for cfg=1:Np
%             figure(43)
%             subplot(1,2,1)
%             semilogy(snrdr,(BBu1(cfg,:)),calcium2{cfg});hold on;
%             subplot(1,2,2)
%             semilogy(snrdr,(BBv1(cfg,:)),calcium2{cfg});hold on;
%             figure(44)
%             semilogy(snrdr,(BB1(cfg,:)),calcium2{cfg});hold on;title('Error Area');xlabel('SNR(dB)');ylabel('Error','Interpreter','Latex');
%             
%         end
%     end
%% ZZB from Pe lower bound 3
% fprintf('\n ZZB3:');
% figure(44)
% for cfg=1:Np
%     [ZZB1(cfg,:),ZZBu1(cfg,:),ZZBv1(cfg,:)]=zzb_intel3(snrdr,txf(:,cfg),tyf(:,cfg),N,ha,vva,A);fprintf(' %d',cfg);
%     if cfg==3
% [ZZB1(cfg,:),ZZBu1(cfg,:),ZZBv1(cfg,:)]=zzb_intel3(snrdr,txf2,tyf2,N,ha,vva,A);
%     end
% end
% figure(43)
% for cfg=1:Np
%     subplot(2,1,1)
%     semilogy(snrdr,5*log10(ZZBu1(cfg,:)),calcium3{cfg});
%     subplot(2,1,2)
%     semilogy(snrdr,5*log10(ZZBv1(cfg,:)),calcium3{cfg});
% end
 %% ZZB from Pe numerical bound 3
% ZZBu2=ZZBu1;ZZBv2=ZZBv1;
fprintf('\n ZZB3b:');
figure(40)
for cfg=1:Np
    subplot(1,3,cfg)
        if cfg<3

    [ZZBi2(cfg,:),ZZBu2(cfg,:),ZZBv2(cfg,:),c_ZZB]=zzb_intel7(snrdr,txf(:,cfg),tyf(:,cfg),N,ha,vva,A,Ns,pb_case);fprintf(' %d',cfg);
        else
            [ZZBi2(cfg,:),ZZBu2(cfg,:),ZZBv2(cfg,:),c_ZZB]=zzb_intel7(snrdr-10*log10(length(txf2)/length(txf(:,1))),txf2,tyf2,N,ha,vva,A,Ns,pb_case);fprintf(' %d',cfg);
    end
end
for cfg=1:Np
    figure(43)
    subplot(2,1,1)
    plot(snrdr,10*log10(ZZBu2(cfg,:)),calcium2{cfg});
    subplot(2,1,2)
    hold on;plot(snrdr,10*log10(ZZBv2(cfg,:)),calcium2{cfg});%hold on;
     figure(44)
    plot(snrdr,10*log10(ZZBi2(cfg,:).^2),calcium2{cfg});hold on;title('Error Area');xlabel('SNR(dB)');ylabel('Error','Interpreter','Latex');
figure(45)    
semilogy(snrdr,5*log10(ZZBu2(cfg,:).^2+ZZBv2(cfg,:).^2),calcium2{cfg});title('Overall Error');xlabel('SNR(dB)');ylabel('Error','Interpreter','Latex');
end
% return;
 %% Numerical MSE
 fprintf('\n MLE:');
 figure(40)
for cfg=1:Np
    subplot(1,3,cfg)
        if cfg<3

    [MLEu(cfg,:),MLEv(cfg,:),MLEi(cfg,:),c_MLE]=DML_MSE(snrdr,txf(:,cfg),tyf(:,cfg),N,Ns,0,pb_case);fprintf(' %d',cfg);
    title(config(cfg));grid on;
        else
            [MLEu(cfg,:),MLEv(cfg,:),MLEi(cfg,:),c_MLE]=DML_MSE(snrdr-10*log10(length(txf2)/length(txf(:,1))),txf2.',tyf2.',N,Ns,0,pb_case);fprintf(' %d',cfg);
    end
end
% legend([c_CRB(1),c_ZZB(1),c_MLE(1)],'CRB','ZZB','MLE','location','northeast');

for cfg=1:Np
    figure(43)
    subplot(2,1,1)
    plot(snrdr,10*log10(MLEu(cfg,:)),calcium4{cfg});hold on;grid on;
    subplot(2,1,2)
    semilogy(snrdr,10*log10(MLEv(cfg,:)),calcium4{cfg});hold on;grid on;
 figure(44)
    plot(snrdr,10*log10(MLEi(cfg,:)),calcium4{cfg});hold on;grid on;
    title('Error Area');xlabel('SNR(dB)');ylabel('Error','Interpreter','Latex');
     figure(45)
    plot(snrdr,5*log10(MLEu(cfg,:).^2+MLEv(cfg,:).^2),calcium4{cfg});hold on;grid on;
    title('Overall Error');xlabel('SNR(dB)');ylabel('Error','Interpreter','Latex');
end
%%
 figure(43)
legend([strcat(config,' -CRB'),strcat(config,' -ZZB'),strcat(config,' -MLE')],'location','northeast');hold on;
 figure(44)
legend([strcat(config,' -CRB'),strcat(config,' -ZZB'),strcat(config,' -MLE')],'location','northeast');hold on;
 figure(45)
legend([strcat(config,' -CRB'),strcat(config,' -ZZB'),strcat(config,' -MLE')],'location','northeast');hold on;
%%
% calcium={'r--','b-.','g:','k-'};calcium2={'r--x','b-.o','g:+','k-s'};
%% Analytically evaluated ZZB
% figure(53)
% for cfg=1:Np
%     subplot(1,2,1)
%     semilogy(snrdr,CRBu(cfg,:),calcium{cfg});hold on;title('RMSE U');xlabel('SNR(dB)');ylabel('RMSE','Interpreter','Latex');axis([sb,se,1e-4,1]);
%     subplot(1,2,2)
%     semilogy(snrdr,CRBv(cfg,:),calcium{cfg});hold on;title('RMSE V');xlabel('SNR(dB)');ylabel('RMSE','Interpreter','Latex');axis([sb,se,1e-4,1]);
% end
% for cfg=1:Np
%     subplot(1,2,1)
%     semilogy(snrdr,sqrt(ZZBu1(cfg,:)),calcium2{cfg});hold on;
%     subplot(1,2,2)
%     semilogy(snrdr,sqrt(ZZBv1(cfg,:)),calcium2{cfg});hold on;
% end
% legend([strcat(config,' -CRB'),strcat(config,' -ZZB1')],'location','northeast');hold on;
% figure(63)
% for cfg=1:Np
%     subplot(1,2,1)
%     semilogy(snrdr,CRBu(cfg,:),calcium{cfg});hold on;title('RMSE U');xlabel('SNR(dB)');ylabel('RMSE','Interpreter','Latex');axis([sb,se,1e-4,1]);
%     subplot(1,2,2)
%     semilogy(snrdr,CRBv(cfg,:),calcium{cfg});hold on;title('RMSE V');xlabel('SNR(dB)');ylabel('RMSE','Interpreter','Latex');axis([sb,se,1e-4,1]);
% end
% for cfg=1:Np
%     subplot(1,2,1)
%     semilogy(snrdr,sqrt(ZZBu2(cfg,:)),calcium2{cfg});hold on;
%     subplot(1,2,2)
%     semilogy(snrdr,sqrt(ZZBv2(cfg,:)),calcium2{cfg});hold on;
% end
% legend([strcat(config,' -CRB'),strcat(config,' -ZZB2')],'location','northeast');hold on;
%% MLE position error plots
if 0
%Ns=64;%N=10;
figure(73)
for cfg=1:Np
    subplot(1,2,1)
    semilogy(snrdr,CRBu(cfg,:),colr(cfg));hold on;title('RMSE U');xlabel('SNR(dB)');ylabel('RMSE','Interpreter','Latex');axis([sb,se,1e-4,1]);
    subplot(1,2,2)
    semilogy(snrdr,CRBv(cfg,:),colr(cfg));hold on;title('RMSE V');xlabel('SNR(dB)');ylabel('RMSE','Interpreter','Latex');axis([sb,se,1e-4,1]);
end
parfor cfg=1:Np
    [MLEu(cfg,:),MLEv(cfg,:),MLEi(cfg,:)]=DML_MSE(snrdr,txf(:,cfg),tyf(:,cfg),N,Ns,0,pb_case);fprintf(' %d',cfg);
end
for cfg=1:Np
    subplot(1,2,1)
    semilogy(snrdr,sqrt(MLEu(cfg,:)),'--','color',colr(cfg));hold on;
    subplot(1,2,2)
    semilogy(snrdr,sqrt(MLEv(cfg,:)),'--','color',colr(cfg));hold on;
end
% parfor cfg=1:Np
%     [MLEu(cfg,:),MLEv(cfg,:),MLEi(cfg,:)]=DML_MSE(snrdr,txf(:,cfg),tyf(:,cfg),N,Ns,1,pb_case);fprintf(' %d',cfg);
% end
% for cfg=1:Np
%         subplot(1,2,1)
%     semilogy(snrdr,sqrt(MLEu(cfg,:)),'-.','color',colr(cfg));hold on;
%     subplot(1,2,2)
%     semilogy(snrdr,sqrt(MLEv(cfg,:)),'-.','color',colr(cfg));hold on;
% end
% parfor cfg=1:Np
%     [MLEu(cfg,:),MLEv(cfg,:),MLEi(cfg,:)]=DML_MSE(snrdr,txf(:,cfg),tyf(:,cfg),N,Ns,2,pb_case);fprintf(' %d',cfg);
% end
% for cfg=1:Np
%     subplot(1,2,1)
%     semilogy(snrdr,sqrt(MLEu(cfg,:)),':','color',colr(cfg));hold on;
%     subplot(1,2,2)
%     semilogy(snrdr,sqrt(MLEv(cfg,:)),':','color',colr(cfg));hold on;
% end
legend([strcat(config,' -CRB'),strcat(config,' -MLE'),strcat(config,' -MLE(\sigma=1mm)'),strcat(config,' -MLE(\sigma=2mm)')],'location','northeast');hold on;
end
%%
U = [-0.5:5e-3:0.5] ;V = U;
% [SU2,SV2] = meshgrid(V,U);
SU2=SU;SV2=SV;
for cfg=1:Np
    figure(7)
    subplot(3,2,2*(cfg-1)+1)
    hold on;grid on;stdmm=0;
    if cfg<3
        M=size(txf,1);
        scatter(txf(:,cfg)-mean(txf(:,cfg)),tyf(:,cfg)-mean(tyf(:,cfg))+randn(M,1)*stdmm/5/sqrt(2),'ks');
        pattern_mult(:,:,cfg)=bpmra_lite(txf(:,cfg),tyf(:,cfg),SU2,SV2,zeros(M,1),[]);
        
    else
        M=length(txf2);
        scatter(txf2-mean(txf2),tyf2-mean(tyf2),'ks','SizeData',10);
        pattern_mult(:,:,cfg)=bpmra_lite(txf2,tyf2,SU2,SV2,zeros(M,1),[]);
    end
    xlabel('x - \lambda units');ylabel('y - \lambda units');axis([-8, 8, -8, 10]);
    %     text(-8,8,config(cfg),'FontSize',14);
    title(config(cfg))
    %% Beampattern
    pattern_multa(:,:,cfg)=mag2db(abs(pattern_mult(:,:,cfg)))+pbi*pb_case;% Add Scan Loss
    %         GD(cfg)=gd(db2mag(pattern_multa(:,:,cfg)),SU,SV);pattern_multa(:,:,cfg)=pattern_multa(:,:,cfg)-max(max(pattern_multa(:,:,cfg)))+GD(cfg);
    pattern_multa(:,:,cfg)=pattern_multa(:,:,cfg)-max(max(pattern_multa(:,:,cfg)))+26;
    figure(7)
    subplot(3,2,2*cfg)
    surf(SU2,SV2,max(pattern_multa(:,:,cfg),0),'linestyle','none');view(2);colorbar;colormap(flipud(colormap));
    hold on;
    %             title([config(cfg),sprintf(' G_D=%2.2f dBi',GD(cfg))]);
    xlabel('U');ylabel('V');zlabel('Gain dB');
    [cost0,msll0,beamw0,eccen0,GD0]= costs_lite1000(pattern_multa(:,:,cfg),mask2,SU2,SV2,pbi);
    %     text(0.2,-0.6,30,sprintf('MSLL:%0.1f dB\n BW:%0.1f deg\n ecc: %0.2f',msll0,beamw0,eccen0),'color','r');
    title(sprintf('MSLL:%0.1f dB, BW:%0.1f deg, ecc: %0.2f',msll0,beamw0,eccen0));
end
%%
mkdir(['Results2/',date_str]);
saveas(figure(40),['Results2/',date_str,'/fig_error_ellipse.fig']);
saveas(figure(43),['Results2/',date_str,'/fig_min_max.fig']);
saveas(figure(44),['Results2/',date_str,'/fig_overall.fig']);
saveas(figure(45),['Results2/',date_str,'/fig_area.fig']);
saveas(figure(7),['Results2/',date_str,'/fig_configBP.fig']);
%%
% figure(8)
% posa=imp_samples(Ns,0.25);
% scatter(posa(1,:),posa(2,:),[],posa(3,:));colorbar;