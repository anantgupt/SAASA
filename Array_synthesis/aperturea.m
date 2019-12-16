function area= aperturea(tx,ty,Nres)
    ang=linspace(0,pi*2*(1-1/Nres),Nres);
    txf=[tx+0.75;tx+0.75;tx-0.75;tx-0.75];tyf=[ty+0.9;ty-0.9;ty+0.9;ty-0.9];
    for i=1:Nres
        apt=(txf+1j*tyf)*exp(-1j*ang(i));
    ap(i)=max(real(apt))-min(real(apt));
    end
%     figure(1)
%     clf;
    fill(ap/2.*cos(ang),ap/2.*sin(ang),'y');axis([-15,15,-15,15]);grid on;hold on;
%     scatter(txf-mean(txf),tyf-mean(tyf),'rx');
    area=trapz((ap./2).^2)*(ang(2)-ang(1))/2;
end