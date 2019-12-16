function [cost,costi]=best_subtree_combs3(pattern,Ni,p,tx,ty)
global subp;global rho;global SU;global SV;global mask2;global flag;
global txc;global tyc;global ncid;global nc;global Np;global Nends;global h;global minslla;global beamwa;global eccena;global gda;
global tx_data;global ty_data;
% fprintf('%d',p);
Nbranches=nc{p}(Ni);
Nstart=ncid{p}(Ni);
p=p+1;
if Nbranches>0
    costa=-inf*ones(Nbranches,1);costai=-inf*ones(Nbranches,1);mslla=costa;bwa=mslla;eca=mslla;
    txa=txc{p}(Nstart:Nstart+Nbranches-1);tya=tyc{p}(Nstart:Nstart+Nbranches-1);
    if p<Np
        Nend=Nbranches;
        for i=1:Nbranches
            txtemp=txa(i);
            tytemp=tya(i);
            patternd=pattern+exp(1j*(2*pi*rho*(txtemp*(SU)+tytemp*(SV))));
            [costa(i),costai(i)]=best_subtree_combs3(patternd,i+Nstart-1,p,[tx,txtemp],[ty,tytemp]);
            if flag==1
                Nend=i;
                break;
            end
        end
        [cost,optj]=min(costa(1:Nend));
        costi=costai(optj);
    else
        
        for i=1:Nbranches
            txtemp=txa(i);
            tytemp=tya(i);
            patternd=pattern+exp(1j*(2*pi*rho*(txtemp*(SU)+tytemp*(SV))));
            [costa(i),mslla(i),bwa(i),eca(i),gd(i)]=costs_lite2(mag2db(abs( patternd.*subp )),mask2);
            tx_data(Nstart+i-1,:)=[tx,txtemp];ty_data(Nstart+i-1,:)=[ty,tytemp];
        end
        minslla(Nstart:Nstart+Nbranches-1)=mslla;
        beamwa(Nstart:Nstart+Nbranches-1)=bwa;
        eccena(Nstart:Nstart+Nbranches-1)=eca;
        gda(Nstart:Nstart+Nbranches-1)=gd;
        [cost,optj]=min(costa);
        costi=optj+Nstart-1;
        waitbar(double(Nstart+Nbranches-1)/Nends,h);
        if (Nstart+Nbranches) >Nends % was>=
            flag=1;
        end
%         disp(size(mslla));disp(minsll);disp(minslli);disp(Ni);disp(Nbranches);
    end
else
    cost=inf;
    costi=Nstart;
end

end