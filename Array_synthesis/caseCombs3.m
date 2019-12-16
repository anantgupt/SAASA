%% May 16 (Last day of 26)
% Finds best config (after cyclic & fine optimization) for low MSLL.
%
%%
clear all;close all;
global rho;global SU;global SV;
global subp;global lastel;global mask2;global flag;
global txc; global tyc; global nc;global ncid; global Np;global h;global Nends;global minslla;global beamwa;global eccena;global gda;
global tx_data;global ty_data;global pbi;

Nres=512;Np=8;rho=1.0;Nxr=4;Nyr=4;
beta=0;gamma=0;
img = imread('rfem.png');             %# Load rfem image
set(0,'DefaultFigureWindowStyle','docked');
in_file='combs9A';
out_file=strcat('data_',in_file,'.mat');
%%
yyy=([0:(Nyr-1)]-(Nyr-1)/2)*0.6;xxx=([0:(Nxr-1)]-(Nxr-1)/2)*0.5;
u0=0;v0=0;nbits=-1;
[agrdx,agrdy]=meshgrid(xxx,yyy);
U = [-Nres/2:(Nres/2)-1] ./(Nres/2);V = [-Nres/2:(Nres/2)-1] ./(Nres/2);
[SU,SV] = meshgrid(V,U);
indx=find(SU.^2+SV.^2>=1);
mask2=(SU.^2+SV.^2) < 1-2/Nres;
bd=@(x,y,z)(x<y).*y+x.*(x>=y).*(x<z)+(x>=z).*z;%bound between y & z
lastel=@(x)x(end);
msy=['x','.','.','.','.','.','.','o'];
%% Patch beampattern
for u=1:size(SU,1)
    for v=1:size(SU,2)
        pb(v,u)=patchbeam(SU(v,u),SV(v,u));
    end
end
pbi=20*log10(pb);
%% Subarray beampattern
subp=bpmra_lite(agrdx(:),agrdy(:),SU,SV,ones(numel(agrdx),1),indx,rho);
txo=zeros(Np,length(xxx)/2*length(yyy)/2);tyo=txo;roto=txo;mino=zeros(1,length(xxx)/2*length(yyy)/2);
%%
load([in_file,'.mat']);
txc=Q.txc;
tyc=Q.tyc;
nc=Q.nc;
% rotc3=Q.nc;% For chip-size ignorant case
rotc3=Q.rotc3;
%%
for i=1:Np-1
    nctemp=cumsum(nc{i})+1;
%     nctemp(end)=[];    ncid{i}=[1 nctemp];
    nctemp(2:end)=nctemp(1:end-1);nctemp(1)=1;    ncid{i}=nctemp;
end
    pattern=ones(size(SU));Nstart=1;p=1;
    N3c=1;Nends=length(txc{Np});
    h=waitbar(N3c/Nends,'Running...');
    minsll_aray=zeros(size(txc{1}));minslli_aray=minsll_aray;
    txc0=txc{1};tyc0=tyc{1};flag=0;Nend=length(txc{1});
    for i=1:length(txc{1})
        pattern=exp(1j*(2*pi*rho*(txc0(i)*(SU)+tyc0(i)*(SV))));
            [minsll_aray(i),minslli_aray(i)]=best_subtree_combs3(pattern,i,p,txc0(i),tyc0(i));
            if flag==1
                Nend=i;
                break;
            end
    end
    delete(h);
    [minsll,indx_temp]=min(minsll_aray(1:Nend));minslli=minslli_aray(indx_temp);
    %%
        indx=minslli;
            txp(Np)=txc{Np}(indx);typ(Np)=tyc{Np}(indx);rotp(Np)=rem(idivide(rotc3(minslli),3^(Np-1)),3);
            j=Np-1;
        while j>0
            indx=find(indx<=cumsum(nc{j}),1);
            txp(j)=txc{j}(indx);typ(j)=tyc{j}(indx);
            rotp(j)=rem(idivide(rotc3(minslli),3^(j-1)),3);% Number packed rot
            j=j-1;
        end
%%
tx=txp.';ty=typ.';rot=(rotp.'==1)*pi;
%%
        tmp=repmat(tx,[1,16])+ones(Np,1)*agrdx(:).';
    txf=tmp(:);
    tmp=repmat(ty,[1,16])+ones(Np,1)*agrdy(:).';
    tyf=tmp(:);

%%
%     close(v1);
clf(figure(11));
clf(figure(32));
figure(11)
% aperturea(tx,ty,Nres);hold on;
scatter(txf,tyf,'bs');hold on;grid on;
plot(tx,ty,'rx-');title('Array element Positions');legend('element','subarray center','location','southeast');
% for i=1:Np
%     rectangle('Position',[tx(i),ty(i)+rot(i)/pi*3.68-1.84,2,5.08],'facecolor','g');
% end
figure(32)
for i=1:Np
    xw=abs(2*cos(rot(i)))+abs(5.04*sin(rot(i)));
    yw=abs(2*sin(rot(i)))+abs(5.04*cos(rot(i)));
    image([tx(i)+1.32*sin(rot(i))-xw/2 tx(i)+1.32*sin(rot(i))+xw/2],[ty(i)-1.32*cos(rot(i))+yw/2 ty(i)-1.32*cos(rot(i))-yw/2],imcomplement(imrotate(imcomplement(img),rot(i)*180/pi))); hold on;
    axis([-12, 12, -12, 15]);
end
scatter(tx(:),ty(:),'r.')
set(gca,'YDir','normal');xlabel('x - \lambda units');ylabel('y - \lambda units');
%%
rho=1.5;
utrue=0.0;vtrue=0.0;
figure(13)
subplot(1,2,1)
surf(SU,SV,max(0,mag2db(abs( bpmra(txf,tyf,SU,SV,u0,v0,nbits,indx,rho) ))),'linestyle','none');view(2);
hold on; circt=1/rho*exp(1j*(0:0.1:2*pi))-(utrue+vtrue*1j)/rho;plot3(real(circt),imag(circt),44*ones(size(circt)),'r.-');colorbar;
title('Modified Beampattern(MBP)');xlabel('U`');ylabel('V`')
subplot(1,2,2)
surf(SU,SV,max(0,pbi+mag2db(abs( bpmra(txf,tyf,SU,SV,utrue,vtrue,nbits,indx) ))),'linestyle','none');view(2);colorbar;
title(sprintf('Original beampattern steered to %0.1f,%0.1f',utrue,vtrue));
msll_val_final=(costs_lite(pbi+mag2db(abs( bpmra(txf,tyf,SU,SV,utrue,vtrue,nbits,indx) )),mask2));
fprintf('MSLL: %2.2f',msll_val_final);
%%
figure(15)
scatter(Q.tra,Q.deta,[],minslla,'.');colorbar;hold on;
scatter(Q.tra(minslli),Q.deta(minslli),[],'r','x');
xlabel('\lambda_2');ylabel('\lambda_1');text(Q.tra(minslli),real(Q.deta(minslli)),num2str(msll_val_final),'Color','red','FontSize',14);
%%
figure(16)
scatter3(Q.tra,Q.deta,minslla,[],minslla,'.');colorbar;hold on;
scatter3(Q.tra(minslli),Q.deta(minslli),minsll,[],'r','x');
xlabel('\lambda_2');ylabel('\lambda_1');text(Q.tra(minslli),real(Q.deta(minslli)),minslla(minslli),num2str(msll_val_final),'Color','red','FontSize',14);
%%
pattern_final=pbi+mag2db(abs( bpmra(txf,tyf,SU,SV,utrue,vtrue,nbits,indx) ));
GD=gd(db2mag(pattern_final),SU,SV);
[msll,asll,bwa,msls,ptge,umax,vmax,id,id2,gd2,pm,pm2,ecc]= costs(squeeze(pattern_final),SU,SV,u0,v0,GD);
fprintf('\n G_D=%0.1f MSLL:%0.1f dBi, BW:%0.1f deg^2, ASLL:%0.1f dBi, Ecc: %0.2f\n',GD,msll,bwa,asll,ecc);

%%
figure(19)
% scatter(minslla,eccena,'.');ylabel('Scale Factor, s');xlabel('MSLL')
hist(eccena,500);
%%
Q.minslla=minslla;Q.beamwa=beamwa;Q.eccena=eccena;Q.gda=gda;
Q.tx_data=tx_data;Q.ty_data=ty_data;Q.res=Nres;Q.rot_data=rotp;
matfile(out_file);
save(out_file,'Q');
%%
% load('May24_2018b.mat');
% pp=1;name='Triangle_minbw1';
% P.txa(:,pp)=tx;P.tya(:,pp)=ty;P.rot(:,pp)=rot;P.name=name;
% save May24_2018b.mat P;
%%
% t_data=[tx_data,ty_data];tbl=array2table(t_data);
% gprMdl = fitrgp(tbl,beamwa,'KernelFunction','ardsquaredexponential',...
%       'FitMethod','sr','PredictMethod','fic','Standardize',1);
  %% Show model fit quality
%   ypred = resubPredict(gprMdl);
%   figure();
% plot(tbl.NoShellRings,'r.');
% hold on
% plot(ypred,'b');
% xlabel('x');
% ylabel('y');
% legend({'data','predictions'},'Location','Best');
% axis([0 4300 0 30]);
% hold off;
% L = resubLoss(gprMdl)%regression loss on the training data (resubstitution loss) 
%% 
figure(17)
scatter3(Q.tra,Q.deta,ov_cost(beamwa,minslla,eccena,gda),[],ov_cost(beamwa,minslla,eccena,gda),'.');hold on;
scatter3(Q.tra(minslli),Q.deta(minslli),minsll,[],'r','x');
% textbox(Q.tra(minslli),real(Q.deta(minslli)),minsll,'String',num2str(minsll));
xlabel('\lambda_2');ylabel('\lambda_1');title('Cost Function');
%%
figure(18)
scatter(Q.beamwa,Q.minslla,'.');hold on;xlabel('BW');ylabel('MSLL')
scatter(Q.beamwa(minslli),Q.minslla(minslli),'rx')
