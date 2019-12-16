function [txf2,tyf2]= fill_aperture(txff,tyff,dolx,doly)
txff=txff-mean(txff);
tyff=tyff-mean(tyff);
indc=convhull(txff,tyff);
txc=txff(indc)/dolx;tyc=tyff(indc)/doly;
%interpolate convhull
txci=[];tyci=[];
for i=1:length(indc)-1
    M=length(txci);
    N=2*ceil(max(abs(txc(i+1)-txc(i)),abs(tyc(i+1)-tyc(i))));
        xi=linspace(txc(i),txc(i+1),N+1);
        yi=linspace(tyc(i),tyc(i+1),N+1);
        K=unique(round([xi(1:N).',yi(1:N).']),'rows','stable');
        Nn=size(K,1);
        txci(M+1:M+Nn)=K(:,1);
        tyci(M+1:M+Nn)=K(:,2);
end
% Find y min,max across x-grid
x_max=txci;
y_max=tyci;
i_start=min(x_max)-1;
x_grid=min(x_max):max(x_max);
y_minmax=[];
yl=zeros(size(x_grid));
for i=1:length(x_max)
    indx=x_max(i)-i_start;
    y_minmax(indx,yl(indx)+1)=y_max(i);
    yl(indx)=yl(indx)+1;
end
% Fill y positions in x-grid
txf2=[];tyf2=[];
for i=1:length(x_grid)
    y_2(1)=min(y_minmax(i,1:yl(i)));y_2(2)=max(y_minmax(i,1:yl(i)));
    Mx=length(txf2);
    txf2(Mx+1:Mx+y_2(2)-y_2(1)+1)=x_grid(i)*ones(y_2(2)-y_2(1)+1,1);
    tyf2(Mx+1:Mx+y_2(2)-y_2(1)+1)=y_2(1):y_2(2);
end

txf2=txf2*dolx;
tyf2=tyf2*doly;
% plotting
% figure(23)
% plot(txci*dolx,tyci*doly,'rs-');hold on;
% plot(txc*dolx,tyc*doly,'m*-');grid on;
% scatter(txf2,tyf2,'k.');
% scatter(txff,tyff,'gs');
end