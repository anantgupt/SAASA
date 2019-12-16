function out = ov_cost(beamw,msll,eccen,gd)
alpha=0.5;beta=0.1;gamma=0.0;psi=1;chi=0.5;%20
bw2=beamw.*sqrt(1-eccen.^2);def=1;
switch def
    case 1% Max Beamwidth
        BW=beamw;
    case 2% Area
        BW=beamw.*bw2;
    case 3% MSE BW
        BW=beamw+bw2;
end
out=alpha*BW/50+beta*msll/10+psi*(msll>-10)+psi*(gd<24)+gamma*real(eccen)+chi*gd/27;
end