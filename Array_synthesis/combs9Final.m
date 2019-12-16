%% Combs Non-parallel
% Unique pair of eigenvalues
% May 15, 2018
% Use trace of correlation as extra parameter to get more choices at
% Considering Rotation 2 as free variable
%%
clear all;close all;
Nxr=4;Nyr=4;
dolx=0.2;doly=0.2;
Np=12;rho=1.0;N_init=16+1;
wi=2;hi=5.04;his=2.4;
D_e = sqrt(dolx^2+doly^2);

R_fac1=4/(2*D_e);R_fac2=8/(2*D_e);% was 20
X_max=10;Y_max=10; % was 8
apx0=round(-X_max:dolx:X_max,3);apy0=round(-Y_max:doly:Y_max,3);%Aperture size(To be fixed according to ZZB)
[apy,apx]=meshgrid(apy0,apx0);
%%
yyy=([0:(Nyr-1)]-(Nyr-1)/2)*doly;xxx=([0:(Nxr-1)]-(Nxr-1)/2)*dolx;
u0=0;v0=0;nbits=-1;
[agrdx,agrdy]=meshgrid(xxx,yyy);
tra=[];deta=[];
%%
nc=cell(Np,1);txc=nc;tyc=nc;nc_cum=cell(Np,1);
p=1;txc{1}=single(linspace(min(apx0),0,N_init));tyc{1}=single(linspace(min(apx0),0,N_init));Nlast=N_init;rotc3=uint16(2*ones(1,N_init));%0-fix 0;1-fix pi;2-both allowed
maskx_plus=int16([1:numel(apx)]);apxf=single(apx(maskx_plus));apyf=single(apy(maskx_plus)); % Regular Grid
ap=rand(2,numel(apx));
% apxf=single((ap(1,:)-0.5)*2*X_max);apyf=single((ap(2,:)-0.5)*2*Y_max); % Random Grid
indxrem=logical(maskx_plus);
txmo=single(txc{1});tymo=single(tyc{1});
h=waitbar(p/Np,'Running...');flag=0;tic;
vara=zeros(N_init,1);pdistmo=zeros(N_init,1);% 2nd order and 1st order moment of superarray pdist

Nraw_branches = [N_init, zeros(1, Np-1)]; % Raw branches in iterations
N_prune1 = [N_init, zeros(1, Np-1)]; % Brances left after removing high det
N_prune2 = [N_init, zeros(1, Np-1)]; % Branches left after 1st order Pruning (among siblings)
N_prune3 = [N_init, zeros(1, Np-1)]; % Branches left after 2nd order Pruning (all nodes at same level)
while p<Np
    waitbar(p/Np,h);
    txtemp=single([]);tytemp=single([]);rottemp=uint16([]);
    txci=txc{p};tyci=tyc{p};Next_indx=1;
    deta=[];tra=[];
    nctemp=zeros(Nlast,1,'uint32');
    p=p+1;
    txtree=zeros(p-1,1);tytree=txtree;rotree=uint16(txtree);
    for i=1:Nlast
        %% Domain for current placement
        indx=i;
        ro=rem(idivide(rotc3(indx),3^(p-2),'floor'),3);
        if ro==2
            indxrem1=indxrem & (abs(apxf-txci(indx))>=2 | (apyf-tyci(indx))>=his | (apyf-tyci(indx))<=-hi);
            indxrem2=indxrem & (abs(apxf-txci(indx))>=2 | (apyf-tyci(indx))>=hi | (apyf-tyci(indx))<=-his);
        else
            indxrem1=indxrem & (abs(apxf-txci(indx))>=2 | abs(apyf+1.32-(tyci(indx)+(2*single(ro)-1)*1.32))>=hi);%for pi
            indxrem2=indxrem & (abs(apxf-txci(indx))>=2 | abs(apyf-1.32-(tyci(indx)+(2*single(ro)-1)*1.32))>=hi);%for 0
        end
        txtree(p-1)=txc{p-1}(indx);tytree(p-1)=tyc{p-1}(indx);rotree(p-1)=ro;
        j=p-2;
        while j>0
            indx=find(indx<=nc_cum{j},1);
            ro=rem(idivide(rotc3(i),3^(j-1),'floor'),3);
            if ro==2
                indxrem1=indxrem1 & (abs(apxf-txc{j}(indx))>=2 | (apyf-tyc{j}(indx))>=his | (apyf-tyc{j}(indx))<=-hi);% index, rev placement
                indxrem2=indxrem2 & (abs(apxf-txc{j}(indx))>=2 | (apyf-tyc{j}(indx))>=hi | (apyf-tyc{j}(indx))<=-his);% index, fwd
            else
                indxrem1=indxrem1 & (abs(apxf-txc{j}(indx))>=2 | abs(apyf+1.32-(tyc{j}(indx)+(2*single(ro)-1)*1.32))>=hi);%for pi
                indxrem2=indxrem2 & (abs(apxf-txc{j}(indx))>=2 | abs(apyf-1.32-(tyc{j}(indx)+(2*single(ro)-1)*1.32))>=hi);%for 0
            end
            txtree(j)=txc{j}(indx);tytree(j)=tyc{j}(indx);rotree(j)=ro;%Store posn, rot of subtree
            j=j-1;
        end
        %% rotation correction
        rot_cor1=uint16(0*indxrem);rot_cor2=rot_cor1;rota=rot_cor1+rotc3(i);
        for j=1:p-1
            if rotree(j)==2
                % rotation correction of past placements
                rot_cor1=uint16(0*indxrem);rot_cor2=rot_cor1;
                rot_cor1(indxrem1)=uint16(2*(3^(j-1))*(abs(apxf(indxrem1)-txtree(j))<2 & (apyf(indxrem1)-tytree(j))<his & (apyf(indxrem1)-tytree(j))>(his-2*hi))+...
                    (3^(j-1))*(abs(apxf(indxrem1)-txtree(j))<2 & (tytree(j)-apyf(indxrem1))<hi & (tytree(j)-apyf(indxrem1))>-hi));
                rot_cor2(indxrem2)=uint16(2*(3^(j-1))*(abs(apxf(indxrem2)-txtree(j))<2 & (apyf(indxrem2)-tytree(j))<(2*hi-his) & (apyf(indxrem2)-tytree(j))>-his)+...
                    (3^(j-1))*(abs(apxf(indxrem2)-txtree(j))<2 & (tytree(j)-apyf(indxrem2))<hi & (tytree(j)-apyf(indxrem2))>-hi));
                rota=rota-max(rot_cor1,rot_cor2);
            end
            % correction in rot of pth placement
            roto=rem(idivide(rota,3^(j-1),'floor'),3);
            indxrem1=indxrem1 & (abs(apxf-txtree(j))>=2 | (abs(apyf+2*1.32-(tytree(j)+(2*single(roto))*1.32))>=hi & roto<2 ) | roto==2);%for pi
            indxrem2=indxrem2 & (abs(apxf-txtree(j))>=2 | (abs(apyf-(tytree(j)+(2*single(roto))*1.32))>=hi & roto<2 ) | roto==2);%for 0
        end
        rotaray=(uint16(indxrem1 & indxrem2)+uint16(indxrem1))*(3^(p-1))+rota;
        
        %% Branches for current node
        indxremf=indxrem1 | indxrem2;
        txa=apxf(indxremf);tya=apyf(indxremf);rot_pruned=rotaray(indxremf);
        txm=txmo(i)*(p-1)/p+txa/p;tym=tymo(i)*(p-1)/p+tya/p;% Mean of superarray
        pdist_mat=bsxfun(@minus,txtree,txa).^2+bsxfun(@minus,tytree,tya).^2;% Euclidean dist for new center
        pdistm=pdistmo(i)*(p-2)/p+sum(sqrt(pdist_mat),1)*2/p/(p-1);
        c11=(txa-txm).*(txa-txm)/p;c22=(tya-tym).*(tya-tym)/p;c12=(txa-txm).*(tya-tym)/p;
        var0=vara(i)*(p-2)/p+sum(pdist_mat,1)*2/p/(p-1);
        std=real(sqrt(var0-pdistm.^2));
        for j=1:p-1
            c11=(txtree(j)-txm).^2/p+c11;
            c22=(tytree(j)-tym).^2/p+c22;
            c12=(txtree(j)-txm).*(tytree(j)-tym)/p+c12;
        end
        quad_fac=sqrt((c11-c22).^2+4*c12.^2);
        det0=(c11+c22+quad_fac)*0.5;
        tr0=(c11+c22-quad_fac)*0.5;
        % Prune further
        indxallowed1=(c11.*c22-c12.^2)<200 ;%& (std > 0.8*max(std));% Discard high determinant & uniform confg
        temp=find(indxallowed1);
%         [~,ia,~]=unique(uint16([(R_fac1*det0(temp).'),(R_fac1*tr0(temp).'),R_fac1*(std(temp).')/5]),'stable','rows');
        [~,ia,~]=uniquetol(([(det0(temp).'),(tr0(temp).'),(std(temp).')]),1/R_fac1,'Byrows',true,'DataScale',1);
%         [~,ia,~]=uniquetol(([(det0(temp).'),(tr0(temp).')]),1/R_fac1,'Byrows',true,'DataScale',1);
        indxallowed=temp(ia);
        N_Branches=nnz(indxallowed);
        % SAve number of branches
        Nraw_branches(p) = Nraw_branches(p)+length(det0);
        N_prune1(p) = N_prune1(p) + length(temp); % Brances left after removing high det
        N_prune2(p) = N_prune2(p) + N_Branches; % Branches left after 1st order Pruning (among siblings)

        nctemp(i)=N_Branches; % Branches left after 1st order pruning 
        txtemp(Next_indx:Next_indx+N_Branches-1)=txa(indxallowed);
        tytemp(Next_indx:Next_indx+N_Branches-1)=tya(indxallowed);
        rottemp(Next_indx:Next_indx+N_Branches-1)=rot_pruned(indxallowed);
        txmn(Next_indx:Next_indx+N_Branches-1)=txm(indxallowed);
        tymn(Next_indx:Next_indx+N_Branches-1)=tym(indxallowed);
        deta(Next_indx:Next_indx+N_Branches-1)=single(det0(indxallowed));
        tra(Next_indx:Next_indx+N_Branches-1)=single(tr0(indxallowed));
        varn(Next_indx:Next_indx+N_Branches-1)=single(var0(indxallowed));
        pdistmn(Next_indx:Next_indx+N_Branches-1)=single(pdistm(indxallowed));
%%
        Next_indx=Next_indx+N_Branches;
%% Parallel loop code
%         txtt{i}=txa;tytt{i}=tya;
%         c11t{i}=c11;c12t{i}=c12;c22t{i}=c22;
%         dett{i}=det0;trt{i}=tr0;
    end
    varn(Next_indx:end)=[];pdistmn(Next_indx:end)=[];
        %% Prune branches
%     [C,ia,~]=unique(uint16([(R_fac2*deta.'),(R_fac2*tra.'),(R_fac2*real(sqrt(varn-pdistmn.^2).')/5)]),'stable','rows');%Stable to keep in order, round to avoid floating point comparison error
    [C,ia,~]=uniquetol([(deta.'),(tra.'),real(sqrt(varn-pdistmn.^2).')],1/R_fac2,'Byrows',true,'DataScale',1);%Stable to keep in order, round to avoid floating point comparison error
%     [C,ia,~]=uniquetol([(deta.'),(tra.')],1/R_fac2,'Byrows',true,'DataScale',1);%Stable to keep in order, round to avoid floating point comparison error
    ia=ia.';
    ia=sort(ia);
    txmo=txmn(ia);
    tymo=tymn(ia);
    pdistmo=pdistmn(ia);
    vara=varn(ia);
    %%
    txc{p}=txtemp(ia);
    tyc{p}=tytemp(ia);
    rotc3=uint16(rottemp(ia));
    ncc=cumsum(nctemp);
    Nnew=length(ia);
    N_prune3(p) = Nnew; % Branches left after 2nd order Pruning (all nodes at same level)
    if Nlast>1 
        nctemp(1)=1;
        for i=2:Nlast
            ntempt=find(ia>ncc(i-1),1);
            if isempty(ntempt)
                nctemp(i:Nlast)=(Nnew+1)*ones(Nlast-i+1,1);
                break;
            else
                nctemp(i)=ntempt;
            end
        end
        nctemp(Nlast+1)=Nnew+1;
        nc{p-1}=diff(nctemp);
    else
        nc{p-1}=Nnew;
    end
    nc_cum{p-1}=uint32(cumsum(nc{p-1}));
    Nlast=Nnew; 
    fprintf('Configs in Layer %d : %d <- %d <- %d <- %d\n',p,Nlast, N_prune2(p), N_prune1(p), Nraw_branches(p));
end
delete(h);
Time_elapsed=toc;
fprintf('Time Elapsed=%f\n',Time_elapsed);

%%
matfile('combs9A.mat');
Q.txc=txc;Q.tyc=tyc;Q.nc=nc;
Q.rotc3=rotc3;
Q.tra=tra(ia);Q.deta=deta(ia);Q.vara=real(sqrt(vara-pdistmo.^2));
Q.meta=sprintf('dolx=%1.2f doly=%1.2f R_fac1=%d Rfac2=%d X_max=%2.1f Y_max=%2.1f',dolx,doly,R_fac1,R_fac2,X_max,Y_max);
Q.N_branches = [Nraw_branches; N_prune1; N_prune2; N_prune3];
save combs9A.mat Q;
%% Plot Simulated and expected computational complexity
plot(1:Np,(Q.N_branches(1,:).*[1,10*ones(1,Np-1)]),'bx-');hold on;
plot(1:Np,(Q.N_branches(end,:).*[1,10*ones(1,Np-1)]),'r.-');
N_analyt = (8^4) * (X_max/dolx * Y_max/doly)^4 / 2^9;
S_analyt = (8^3) * (X_max/dolx * Y_max/doly)^3 / 2^9;
plot(1:Np, (ones(1,Np)*N_analyt),'bx--');
plot(1:Np, (ones(1,Np)*S_analyt),'r.--');
legend('Simulated $T_n$', 'Simulated $|\mathcal{C}^n|$', 'Upper Bound, $T_n$', 'Upper Bound, $|\mathcal{C}^n|$');
title('Complexity of Prefix Tree Dictionary Search'); 
ylabel('Operations per iteration');
xlabel('Number of Subarrays Placed (n)');
set(gca, 'YScale', 'log');grid on;
hold off;