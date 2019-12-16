function [cost,iopt] = ov_cost_single(beamw,msll,eccen,gd,Lag_wt,BW_def,BPinfo)
% alpha=0.5;beta=0.1;gamma=0.0;chi=0.0;%20
alpha=Lag_wt(1);beta=Lag_wt(2);gamma=Lag_wt(3);chi=Lag_wt(4);
bw2=beamw.*sqrt(1-eccen.^2);
ApArea=0;%was max(Lag_wt);
switch BW_def
    case 1% Max Beamwidth
        BW=beamw;
    case 2% Area
        BW=beamw.*bw2;
    case 3% MSE BW
        BW=beamw.^2+bw2.^2;%Eigs are sin(\theta)^2
end
%% Normalize objectives
BWn=(BW-BPinfo(1,1))/BPinfo(2,1);
mslln=(msll-BPinfo(1,2))/BPinfo(2,2);
gdn=(gd-BPinfo(1,3))/BPinfo(2,3);
cost=alpha*BWn+beta*mslln+gamma*real(eccen)+chi*gdn;
    
cost=(cost + ApArea*((msll>-10) + (gd<max(gd)-3)));% Add inf if violate constr 
%% Add hierarchical, preference order

end