clear all;close all;
set(0,'DefaultFigureWindowStyle','docked');
%%
in_file='data_combs9A';%combs_dataFinal2b.mat;data_combs9B;combs_data200k;combs_data31;data_combs9A
load([in_file,'.mat']);
img = imread('rfem.png');             %# Load rfem image
% colormap(flipud(colormap));
%% Remove configs with MSLL> -3
% iaccept=find((Q.minslla<-3).*Q.minslla);
iaccept=1:length(Q.minslla);%was -1
N_data=length(iaccept);Nres=400;Np=8;% NOTE: Final Results with Nres=512
yyy=((0:3)-1.5)*0.6;xxx=((0:3)-1.5)*0.5;
nbits=3;u0=0;v0=0;
[agrdx,agrdy]=meshgrid(xxx,yyy);
U = [-Nres/2:(Nres/2)-1] ./(Nres/2);V = [-Nres/2:(Nres/2)-1] ./(Nres/2);
[SU,SV] = meshgrid(V,U);indx=find(SU.^2+SV.^2>=1);
%% Phase optimization parameters
levels = 2^nbits;qlevels = 2.0*pi / levels; % compute quantization levels
NRounds=5;Ntrials=0;plotOrNot=1;msl_ratio=zeros(1,levels);
%%
BW_def=3;% 1: Max BW, 2: BW area, 3:BW MSE
Lag_wt=[0.5,0.5,0.5,0.0];%[BW,MSLL,Ecc,GD]
msll=Q.minslla(iaccept);
eccen=Q.eccena(iaccept);
beamw=Q.beamwa(iaccept);
rotc3=Q.rotc3(iaccept);
gda=Q.gda(iaccept);
% vara=Q.vara(iaccept);
txtest=Q.tx_data(iaccept,:);
tytest=Q.ty_data(iaccept,:);
for id=1:length(iaccept)
    dist0=pdist([txtest(id,:).',tytest(id,:).']);
    vara2(id)=var(dist0);
    mean2(id)=mean(dist0);
    mean3(id)=mean(dist0.^2);
    dist1(id)=var(pdist(txtest(id,:).'));
    dist2(id)=var(pdist(tytest(id,:).'));
    
end
[cost,iopt,BPinfo]=ov_cost_array(beamw,msll,eccen,gda,Lag_wt,BW_def);
cost_opt=cost(iopt);%iopt=1500;

%%
tx=Q.tx_data(iaccept(iopt),:);ty=Q.ty_data(iaccept(iopt),:);
txf=vec(kron(agrdx(:),ones(1,8))+kron(ones(16,1),tx));
tyf=vec(kron(agrdy(:),ones(1,8))+kron(ones(16,1),ty));
rot(Np)=rem(idivide(rotc3(iopt),3^(Np-1)),3);
j=Np-1;
while j>0
    rot(j)=rem(idivide(rotc3(iopt),3^(j-1)),3);% Number packed rot
    j=j-1;
end
rot=(rot.'==1)*pi;
%%
% load(out_file);
% P.txa(:,1)=tx(:);P.tya(:,1)=ty(:);P.rot(:,1)=rot(:);
% Save array
out_file=strcat('array_',in_file,'.mat');
matfile(out_file);
pp=BW_def;name=strcat('Set',in_file(end),'BWdef_',num2str(BW_def));
P.txa(:,pp)=tx(:);P.tya(:,pp)=ty(:);P.rot(:,pp)=rot(:);P.name(:,pp)=name(:);
save(out_file,'P');
%% Patch beampattern
for u=1:size(SU,1)
    for v=1:size(SU,2)
        pb(v,u)=patchbeam(SU(v,u),SV(v,u));
    end
end
pbi=20*log10(pb);
mask2=(SU.^2+SV.^2) < 1-2/Nres;
%% Plotting results
rho=1.5;
utrue=0.0;vtrue=0.0;
figure(14)
% subplot(1,2,2)
pattern_final=mag2db(abs( bpmra(txf,tyf,SU,SV,u0,v0,nbits,indx) ));
GD=gd(db2mag(pbi+pattern_final),SU,SV);
surf(SU,SV,max(0,pattern_final-max(max(pattern_final))+GD),'linestyle','none');view(2);colorbar;
title(sprintf('Original beampattern steered to %0.1f,%0.1f',utrue,vtrue));xlabel('U');ylabel('V')
[msll_val_final,msll_indx]=(costs_lite(pattern_final,mask2));
fprintf('MSLL: %2.2f, Directivity:%0.1f',msll_val_final,GD);
pattern_mbp=mag2db(abs( bpmra(txf,tyf,SU,SV,u0,v0,nbits,indx,rho) ));% use to plot MPB

figure(32)
for i=1:Np
    xw=abs(2*cos(rot(i)))+abs(5.04*sin(rot(i)));
    yw=abs(2*sin(rot(i)))+abs(5.04*cos(rot(i)));
    image([tx(i)+1.32*sin(rot(i))-xw/2 tx(i)+1.32*sin(rot(i))+xw/2],[ty(i)-1.32*cos(rot(i))+yw/2 ty(i)-1.32*cos(rot(i))-yw/2],imcomplement(imrotate(imcomplement(img),rot(i)*180/pi))); hold on;
    axis([-12, 12, -12, 15]);
end
scatter(tx(:),ty(:),'r.')
set(gca,'YDir','normal');xlabel('x - \lambda units');ylabel('y - \lambda units');
[msll_val_final,asll,bwa,msls,ptge,umax,vmax,id,id2,gd2,pm,pm2,ecc]= costs(squeeze(pattern_final),SU,SV,u0,v0,GD);
figure(11)
subplot(1,2,1)
% aperturea(-mean(tx)+tx.',-mean(ty)+ty.',Nres);
hold on;
for i=1:Np
    rectangle('position',[tx(i)-1-mean(tx),ty(i)-mean(ty)+rot(i)/pi*2.68-3.88,2,5.08],'Facecolor',[.97 .97 .97]);
end
scatter(-mean(tx)+txf,-mean(ty)+tyf,'bs');grid on;
plot(-mean(tx)+tx,-mean(ty)+ty,'rx--');title('Subarray Configuration');%legend('element','subarray center','location','southeast');
subplot(1,2,2)
surf(SU,SV,max(0,pattern_final-max(max(pattern_final))+GD),'linestyle','none');view(2);colorbar;
hold on;% circt=1/rho*exp(1j*(0:0.1:6.5))-(utrue+vtrue*1j)/rho;plot3(real(circt),imag(circt),44*ones(size(circt)),'r-');
scatter3(SU(msll_indx),SV(msll_indx),msll_val_final+GD,'rv');
title('Modified Beampattern(MBP)');xlabel('U`');ylabel('V`')
annotation('textbox',[0.02,0.0,0.4,0.05],'String',sprintf('G_D=%0.1f MSLL:%0.1f dBi, ASLL:%0.1f dBi, BW:%0.1f deg^2, Ecc: %0.2f\n',GD,msll_val_final,asll,bwa,ecc));
%% Pareto set figures
txp=tx;typ=ty;
Npert=11;pert_width=0.5;%Perturbation width
subp=bpmra_lite(agrdx(:),agrdy(:),SU,SV,zeros(numel(agrdx),1),indx);
pert_rangex=linspace(-0.5,0.5,Npert);
pert_rangey=linspace(-0.5,0.5,Npert);
[pert_gridx,pert_gridy]=meshgrid(pert_rangex,pert_rangey);
N_SA=5;% #Simulated Annealing rounds
cost_ptb=zeros(Npert);msll_ptb=zeros(Npert);beamw_ptb=zeros(Npert);eccen_ptb=zeros(Npert);GD_ptb=zeros(Npert);
pert_state=NaN(N_SA,Np,5);
%% Video recording
v1 = VideoWriter(['Placement_',in_file,'BW_def',num2str(BW_def),'.mp4'],'MPEG-4');
v1.FrameRate=2;v1.Quality=100;
open(v1);
for sa_rnd=1:N_SA
    howManyElements=1;
    figure(15)%
    clf(figure(15))
    for pert_id=1:Np %subarray to be pertb
        txfp=vec(kron(agrdx(:),ones(1,8))+kron(ones(16,1),txp));
        tyfp=vec(kron(agrdy(:),ones(1,8))+kron(ones(16,1),typ));
        valid=findperts(txp,typ,rot.',xw,yw,pert_id,pert_gridx,pert_gridy);
        
        pattern_mult=bpmra_lite(txfp,tyfp,SU,SV,zeros(size(txfp)),indx);
        pattern_mult_oneless=pattern_mult-subp.*exp(1j*2*pi*(txp(pert_id)*SU+typ(pert_id)*SV));
        tx_pertb=txp(pert_id)+pert_gridx;
        ty_pertb=typ(pert_id)+pert_gridy;
        for pindx1=1:size(pert_gridx,1)
            for pindx2=1:size(pert_gridx,2)
                if isnan(valid(pindx1,pindx2))
                    cost_ptb(pindx1,pindx2)=NaN;
                    msll_ptb(pindx1,pindx2)=NaN;
                    beamw_ptb(pindx1,pindx2)=NaN;
                    eccen_ptb(pindx1,pindx2)=NaN;
                    GD_ptb(pindx1,pindx2)=NaN;
                else
                    pattern_tmp=pattern_mult_oneless+subp.*exp(1j*2*pi*((tx_pertb(pindx1,pindx2))*SU+(ty_pertb(pindx1,pindx2))*SV));
                    
                    [cost_ptb(pindx1,pindx2),msll_ptb(pindx1,pindx2),beamw_ptb(pindx1,pindx2),eccen_ptb(pindx1,pindx2),GD_ptb(pindx1,pindx2)]=...
                        costs_lite1000(mag2db(abs(pattern_tmp)),mask2,SU,SV,pbi,Lag_wt,BW_def,BPinfo);
                end
            end
        end
        inoptb=ceil(numel(tx_pertb)/2);
        [red_cost,red_id]=min(cost_ptb(:));
        txp(pert_id)=tx_pertb(red_id);
        typ(pert_id)=ty_pertb(red_id);
        pert_state(sa_rnd,howManyElements,:) = [red_cost,msll_ptb(red_id),beamw_ptb(red_id),eccen_ptb(red_id),GD_ptb(red_id)];
        howManyElements=howManyElements+1;
        if 0
            figure(16)
            subplot(2,2,1)
            surf(tx_pertb,ty_pertb,log(cost_ptb),'linestyle','none');colorbar;hold on;view(2)
            scatter3(tx_pertb(inoptb),ty_pertb(inoptb),log(cost_ptb(inoptb)),'k','x');title('Cost');
            xlabel('x');ylabel('y');axis([-8,5,-9,6]);
            text(tx_pertb(inoptb),ty_pertb(inoptb),log(max(cost_ptb(:))),num2str(cost_ptb(inoptb)),'Color','red','FontSize',11);
%                 text(-7,-7,log(cost_ptb(inoptb)),num2str(cost_ptb(inoptb)),'Color','red','FontSize',14);
            subplot(2,2,2)
            surf(tx_pertb,ty_pertb,GD_ptb,'linestyle','none');colorbar;hold on;view(2)
            scatter3(tx_pertb(inoptb),ty_pertb(inoptb),GD_ptb(inoptb),'k','x');title('Directivity');
            xlabel('x');ylabel('y');axis([-8,5,-9,6]);
            text(tx_pertb(inoptb),ty_pertb(inoptb),max(GD_ptb(:)),num2str(GD_ptb(inoptb)),'Color','red','FontSize',11);
%             text(-7,-7,GD_ptb(inoptb),num2str(GD_ptb(inoptb)),'Color','red','FontSize',14);
            subplot(2,2,3)
            surf(tx_pertb,ty_pertb,msll_ptb,'linestyle','none');colorbar;hold on;view(2)
            scatter3(tx_pertb(inoptb),ty_pertb(inoptb),msll_ptb(inoptb),'k','x');title('Max Sidelobe Level');
            xlabel('x');ylabel('y');axis([-8,5,-9,6]);
            text(tx_pertb(inoptb),ty_pertb(inoptb),max(msll_ptb(:)),num2str(msll_ptb(inoptb)),'Color','red','FontSize',11);
%             text(-7,-7,msll_ptb(inoptb),num2str(msll_ptb(inoptb)),'Color','red','FontSize',14);
            subplot(2,2,3)
            surf(tx_pertb,ty_pertb,beamw_ptb,'linestyle','none');colorbar;hold on;view(2)
            scatter3(tx_pertb(inoptb),ty_pertb(inoptb),beamw_ptb(inoptb),'k','x');title('BeamWidth');
            xlabel('x');ylabel('y');axis([-8,5,-9,6]);
            text(tx_pertb(inoptb),ty_pertb(inoptb),max(beamw_ptb(:)),num2str(beamw_ptb(inoptb)),'Color','red','FontSize',11);
%             text(-7,-7,beamw_ptb(inoptb),num2str(beamw_ptb(inoptb)),'Color','red','FontSize',14);
        end
        if 1
%             figure(15)
            subplot(3,4,[1 2 5 6])
            for i=1:Np
                xw=abs(2*cos(rot(i)))+abs(5.04*sin(rot(i)));
                yw=abs(2*sin(rot(i)))+abs(5.04*cos(rot(i)));
                image([txp(i)+1.32*sin(rot(i))-xw/2 txp(i)+1.32*sin(rot(i))+xw/2],[typ(i)-1.32*cos(rot(i))+yw/2 typ(i)-1.32*cos(rot(i))-yw/2],imcomplement(imrotate(imcomplement(img),rot(i)*180/pi))); hold on;
                axis tight;grid on;
            end
            scatter(tx(:),ty(:),'r.')
            set(gca,'YDir','normal');title('Array Placement')
            scatter(tx_pertb(~isnan(valid)),ty_pertb(~isnan(valid)),'y.');axis off; 
            subplot(3,4,[3 4 7 8])
            surf(SU,SV,max(10,mag2db(abs(pattern_mult))-max(max(mag2db(abs(pattern_mult))))+GD_ptb(inoptb)),'linestyle','none');view(2);colorbar;
            [msll_val_final,msll_indx]=(costs_lite(mag2db(abs(pattern_mult)),mask2));
            hold on; scatter3(SU(msll_indx),SV(msll_indx),msll_val_final+GD,'rv');axis([-0.5,0.5,-0.5,0.5]);
            title('Beampattern');axis off;
%             annotation('textbox',[0.02,0.0,0.4,0.05],'String',sprintf('G_D=%0.1f MSLL:%0.1f dBi, BW:%0.1f deg^2, Ecc: %0.2f\n',GD_ptb(inoptb)),msll_val_final,beamw_ptb(inoptb),eccen_ptb(inoptb));
subplot(3,4,9)
plot(squeeze(pert_state(:,:,1)).');title('Objective Cost');grid on;axis tight;
subplot(3,4,10)
plot(squeeze(pert_state(:,:,2)).');title('Maximum Sidelobe');grid on;axis tight;
subplot(3,4,11)
plot(squeeze(pert_state(:,:,3)).');title('BeamWidth');grid on;axis tight;
subplot(3,4,12)
plot(squeeze(pert_state(:,:,5)).');title('Directivity');grid on;axis tight;
        end
        F2((sa_rnd-1)*Np+pert_id)=getframe(figure(15));pause(0.01);
    end
    
end
writeVideo(v1,F2);
close(v1);
%%
figure(35)
subplot(2,2,1)
plot(squeeze(pert_state(:,:,1)).');title('Cost');grid on;
subplot(2,2,2)
plot(squeeze(pert_state(:,:,2)).');title('MSLL');grid on;
subplot(2,2,3)
plot(squeeze(pert_state(:,:,3)).');title('BeamWidth');grid on;
subplot(2,2,4)
plot(squeeze(pert_state(:,:,5)).');title('Gd');grid on;
%% Evaluate cost after SA
 txfp=vec(kron(agrdx(:),ones(1,8))+kron(ones(16,1),txp));
tyfp=vec(kron(agrdy(:),ones(1,8))+kron(ones(16,1),typ));
pattern_final=mag2db(abs( bpmra_lite(txfp,tyfp,SU,SV,zeros(size(txfp)),indx) ));
GD=gd(db2mag(pbi+pattern_final),SU,SV);
[msll_val_final,asll,bwa,msls,ptge,umax,vmax,id,id2,gd2,pm,pm2,ecc]= costs(squeeze(pattern_final),SU,SV,u0,v0,GD);
[msll_val_final,msll_indx]=(costs_lite(pattern_final,mask2));
pattern_mbp=mag2db(abs( bpmra(txfp,tyfp,SU,SV,u0,v0,nbits,indx,rho) ));
figure(12)
subplot(1,2,1)
% aperturea(-mean(tx)+tx.',-mean(ty)+ty.',Nres);
hold on;
for i=1:Np
    rectangle('position',[txp(i)-1-mean(txp),typ(i)-mean(typ)+rot(i)/pi*2.68-3.88,2,5.08],'Facecolor',[.97 .97 .97]);
end
scatter(-mean(txp)+txfp,-mean(typ)+tyfp,'bs');grid on;
plot(-mean(txp)+txp,-mean(typ)+typ,'rx--');title('Subarray Configuration');%legend('element','subarray center','location','southeast');
subplot(1,2,2)
surf(SU,SV,max(0,pattern_final-max(max(pattern_final))+GD),'linestyle','none');view(2);colorbar;
hold on; %circt=1/rho*exp(1j*(0:0.1:6.5))-(utrue+vtrue*1j)/rho;plot3(real(circt),imag(circt),44*ones(size(circt)),'r-');
scatter3(SU(msll_indx),SV(msll_indx),msll_val_final+GD,'rv');
title('Modified Beampattern(MBP)');xlabel('U`');ylabel('V`')
annotation('textbox',[0.02,0.0,0.4,0.05],'String',sprintf('G_D=%0.1f MSLL:%0.1f dBi, ASLL:%0.1f dBi, BW:%0.1f deg^2, Ecc: %0.2f\n',GD,msll_val_final,asll,bwa,ecc));
%% Save local refinement
if 1
    load(out_file);
pp=BW_def;name=strcat('Set',in_file(end),'BWdef_',num2str(BW_def));
P.txa(:,pp)=txp(:);P.tya(:,pp)=typ(:);P.rot(:,pp)=rot(:);P.name(:,pp)=name(:);
P.msll(BW_def)=msll_val_final;P.bwa(BW_def)=bwa;
P.Lag_wt(:,pp)=Lag_wt;
save(out_file,'P');
end

if 0
    figure(16)%
    subplot(1,3,1)
    scatter(Q.tra(iaccept),Q.deta(iaccept),[],msll(iaccept),'.');colorbar;hold on;
    scatter(Q.tra(iopt),Q.deta(iopt),[],msll(iopt),'r','x');
    xlabel('\lambda_2');ylabel('\lambda_1');text(Q.tra(iopt),real(Q.deta(iopt)),num2str(msll_val_final),'Color','red','FontSize',14);
    subplot(1,3,2)
    scatter(Q.tra(iaccept),Q.deta(iaccept),[],beamw,'.');colorbar;hold on;
    scatter(Q.tra(iopt),Q.deta(iopt),[],bwa,'r','x');
    xlabel('\lambda_2');ylabel('\lambda_1');text(Q.tra(iopt),real(Q.deta(iopt)),num2str(bwa),'Color','red','FontSize',14);
    subplot(1,3,3)
    scatter(Q.tra(iaccept),Q.deta(iaccept),[],real(eccen),'.');colorbar;hold on;
    scatter(Q.tra(iopt),Q.deta(iopt),[],real(eccen(iopt)),'r','x');
    xlabel('\lambda_2');ylabel('\lambda_1');text(Q.tra(iopt),real(Q.deta(iopt)),num2str(real(eccen(iopt))),'Color','red','FontSize',14);
end
%% Overall Cost function
figure(17)
% subplot(1,2,1)
% scatter3(Q.tra(iaccept),Q.deta(iaccept),cost,[],cost,'b.');hold on;
% scatter3(Q.tra(iaccept(iopt)),Q.deta(iaccept(iopt)),cost_opt,[],'r','d');
% % textbox(Q.tra(minslli),real(Q.deta(minslli)),minsll,'String',num2str(minsll));
% xlabel('\lambda_2');ylabel('\lambda_1');zlabel('f');title(sprintf('Cost Function vs Latent space variabels, BW definition=%d',BW_def));
% subplot(1,2,2)
bw2=beamw.*sqrt(1-eccen.^2);
switch BW_def
    case 1% Max Beamwidth
        BW=beamw;
    case 2% Area
        BW=beamw.*bw2;
    case 3% MSE BW
        BW=beamw.^2+bw2.^2;%Eigs are sin(\theta)^2
end
subplot(1,2,1)
scatter(beamw,msll,'b.');hold on;xlabel('BW_1');ylabel('MSLL')
scatter(beamw(iopt),msll(iopt),'rd');grid on;title(sprintf('Objective Tradeoffs BW definition=%d',BW_def));
subplot(1,2,2)
scatter(BW,msll,'b.');hold on;xlabel('BW_{def}');ylabel('MSLL')
scatter(BW(iopt),msll(iopt),'rd');grid on;title(sprintf('Objective Tradeoffs BW definition=%d',BW_def));

%%
% h(1)=figure(1);
% h(2)=figure(2);
% savefig(h,'grpfit_bw_pdist_expsq_dB.fig');
% close(h);
%% Try sequential phase optimization
starea=find(abs(SU).^2+ abs(SV).^2<0.25);starea=starea(round(rand(1,Ntrials)*numel(starea)));ua=[0.3,0.86,-0.7];va=[0.1,0,-0.5];
for trialNum = 1:Ntrials
    u0=SU(starea(trialNum));%u0=0;u0=ua(trialNum);
    utrue(trialNum)=u0;
    v0=SV(starea(trialNum));%v0=0;v0=va(trialNum);
    vtrue(trialNum)=v0;
    theta0(trialNum)=asin(sqrt(u0^2+v0^2));
    sintheta=txf * u0 +tyf*v0; % Ideal phases
    sinthetaq= round(sintheta* levels ) * qlevels; % Naive rounding
    %% REsults prior to optimization
    pattern=mag2db(abs(bpmra_lite(txf,tyf,SU,SV,sinthetaq,indx)));
    [cost0,msll0,beamw0,eccen0,GD0]= costs_lite1000(pattern,mask2,SU,SV,pbi,Lag_wt,1,BPinfo);% To extract bw1, bw2
    [msll_val_final,msll_indx]=(costs_lite(pattern,mask2));
    fprintf('\n Pre-Phase: G_D=%0.1f dBi, MSLL:%0.1f, BW1:%0.1f deg, BW2=%0.1f Ecc: %0.2f\n',GD0,msll0,beamw0,beamw0.*sqrt(1-eccen0.^2),eccen0);
    figure(13)
    subplot(2,Ntrials,trialNum)
    surf(SU,SV,max(0,pattern-max(pattern(:))+GD0),'linestyle','none');view(2);
    hold on; scatter3(u0,v0,GD0,'r+');colorbar;
    scatter3(SU(msll_indx),SV(msll_indx),msll_val_final+GD0,'rv');
    annotation('textbox',[0.01+(trialNum-1)/3,0.5,0.3,0.03],'String',sprintf('G_D=%0.1f dBi, MSLL:%0.1f dB, Bw1:%0.1f deg, BW2:%0.1f deg, Ecc: %0.2f\n',GD0,msll0,beamw0,beamw0.*sqrt(1-eccen0.^2),eccen0));
    title('Pre-Phase: Modified Beampattern(MBP)');xlabel('U`');ylabel('V`');axis([-1,1,-1,1]);
    %%
    sintheta_seq = sinthetaq;
    tempStorage = zeros(16*Np,5);
    pattern_mult=bpmra_lite(txf,tyf,SU,SV,sinthetaq,indx);
    pattern_mult_seq=pattern_mult;
    fprintf('\n MSLL after ');
    for roundID = 1:NRounds
        howManyElements = 0;
        for elementID = randperm(16*Np)
            howManyElements = howManyElements + 1;
            xtemp=txf(elementID);ytemp=tyf(elementID);
            pattern_tmp=pattern_mult_seq-exp(1j*(2*pi*(xtemp*SU+ytemp*SV)-sintheta_seq(elementID)));
            %% Sequentially optimize cost
            for i=1:levels
                pattern_i=mag2db(abs(pattern_tmp+exp(1j*(2*pi*(xtemp*(SU)+ytemp*(SV))-(i-1)*qlevels))));%NOTE: rem pbi
                [~,maxi]=max(pattern_i(:));
                [msl_ratio(i),msll_i(i),bw_i(i),ecc_i(i),gd_i(i)]=costs_lite1000(pattern_i,mask2,SU,SV,pbi,Lag_wt,BW_def,BPinfo);
                msl_ratio(i)=msl_ratio(i)+max(Lag_wt)*10*sqrt((SU(maxi)-u0)^2+(SV(maxi)-v0)^2);% Steering direction doesn't deviate too much
            end
            [maxmsl,maxIndex] = min(msl_ratio);
            tempStorage(howManyElements,:) = [maxmsl,msll_i(maxIndex),bw_i(maxIndex),ecc_i(maxIndex),gd_i(maxIndex)];
            sintheta_seq(elementID) = (qlevels)*(maxIndex-1);
            pattern_mult_seq=pattern_tmp+exp(1j*(2*pi*(xtemp*(SU)+ytemp*(SV))-(qlevels)*(maxIndex-1)));
        end
        %%
        if plotOrNot
            figure(34)
            for gh=1:5
                subplot(5,Ntrials,trialNum+(gh-1)*Ntrials)
                plot(tempStorage(:,gh));hold on;
                brkpoint = 0;
            end
        end
        fprintf('Rd%d: %2.2f, ',roundID,maxmsl);
    end
    %% Results after optimization
    pattern=mag2db(abs(bpmra_lite(txf,tyf,SU,SV,sintheta_seq,indx)));
    [cost0,msll0,beamw0,eccen0,GD1]= costs_lite1000(pattern,mask2,SU,SV,pbi,Lag_wt,1,BPinfo);
    [msll_val_final,msll_indx]=(costs_lite(pattern,mask2));
    fprintf('\n Post-Phase: G_D=%0.1f MSLL:%0.1f dBi, BW1:%0.1f deg, BW2=%0.1f Ecc: %0.2f\n',GD1,msll0,beamw0,beamw0.*sqrt(1-eccen0.^2),eccen0);
    figure(13)
    subplot(2,Ntrials,trialNum+Ntrials)
    surf(SU,SV,max(0,pattern-max(max(pattern))+GD1),'linestyle','none');view(2);
    hold on; scatter3(u0,v0,GD1,'r+');colorbar;
    scatter3(SU(msll_indx),SV(msll_indx),msll_val_final+GD1,'rv');
    annotation('textbox',[0.01+(trialNum-1)/3,0.0,0.3,0.03],'String',sprintf('G_D=%0.1f dBi, MSLL:%0.1f dB, Bw1:%0.1f deg, BW2:%0.1f deg, Ecc: %0.2f\n',GD1,msll0,beamw0,beamw0.*sqrt(1-eccen0.^2),eccen0));
    title('Post-Phase: Modified Beampattern(MBP)');xlabel('U`');ylabel('V`');axis([-1,1,-1,1]);
end
return;
%% MRA analysis
pdx=pdist(tx.');pdy=pdist(ty.');
figure(41)
subplot(2,2,1);
% scatter(pdx,pdy);
surf(SU,SV,max(0,mag2db(abs(bpmra_lite(tx,ty,SU,SV,zeros(size(tx)),indx)))-0),'linestyle','none');view(2)
tempxsub=[];tempysub=[];
for i=1:length(agrdx(:))
    tempxsub(i,:)=agrdx(i)-agrdx(:);
    tempysub(i,:)=agrdy(i)-agrdy(:);
end
% for i=1:length(tx(:))
% tempx(i,:)=tx(i)-tx(:);
% tempy(i,:)=ty(i)-ty(:);
% end
gdres=0.5;
xedges1=(round(min(tempxsub(:))/gdres)*gdres-gdres/2):gdres:(round(max(tempxsub(:))/gdres)*gdres+gdres/2);
yedges1=(round(min(tempysub(:))/gdres)*gdres-gdres/2):gdres:(round(max(tempysub(:))/gdres)*gdres+gdres/2);

[Nsub,pdsubx,pdsuby]=histcounts2(tempxsub(:),tempysub(:),xedges1,yedges1);
[xgrdsub,ygrdsub]=meshgrid((pdsubx(1:end-1)+pdsubx(2:end))/2,(pdsuby(1:end-1)+pdsuby(2:end))/2);xgrdsub=xgrdsub.';ygrdsub=ygrdsub.';
subplot(2,2,2)
h1=histogram2(tempxsub(:),tempysub(:),xedges1,yedges1);
h1.FaceColor='flat';
% surf(SU,SV,max(0,mag2db(abs(bpmra_lite(txf,tyf,SU,SV,zeros(size(txf)),indx)))-20),'linestyle','none');view(2)
tempx=[];tempy=[];
subplot(2,2,3)
for i=1:length(txf(:))
    tempx(i,:)=txf(i)-txf(:);
    tempy(i,:)=tyf(i)-tyf(:);
end
xedges=(round(min(tempx(:))/gdres)*gdres-gdres/2):gdres:(round(max(tempx(:))/gdres)*gdres+gdres/2);
yedges=(round(min(tempy(:))/gdres)*gdres-gdres/2):gdres:(round(max(tempy(:))/gdres)*gdres+gdres/2);

h=histogram2(tempx(:),tempy(:),xedges,yedges);h.FaceColor='flat';
axis([-11,11,-11,11]);
[N,pdx,pdy]=histcounts2(tempx(:),tempy(:),xedges,yedges);
[xgrd,ygrd]=meshgrid((pdx(1:end-1)+pdx(2:end))/2,(pdy(1:end-1)+pdy(2:end))/2);xgrd=xgrd.';ygrd=ygrd.';
Nf=N;
for i=1:numel(xgrdsub)
    ifound=find((abs((xgrd(:)-xgrdsub(i)))+abs(ygrd(:)-ygrdsub(i))) <gdres/2);
    if numel(ifound)>1
        disp('Errr');
    end
    if ~isempty(ifound)
        Nf(ifound)=Nf(ifound)-Nsub(i)*Np;
        %         fprintf('%d, %d\n',ifound,i);
    end
end
emtyindx=find(N==0);fillindx=find(N~=0);
subplot(2,2,4)
% surf(xgrd,ygrd,log(1+Nf),'linestyle','none');view(2);axis([-11,11,-11,11]);
surf(xgrd,ygrd,Nf,'linestyle','none');view(2);axis([-11,11,-11,11]);
%%
% function out=vec(in)
% out=in(:);
% end