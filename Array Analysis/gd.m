function out = gd(bp,ua,va)
out=0;
R=ua.^2+va.^2;
secc=zeros(size(R));
m=max(max(bp));
for i=1:size(bp,2)
    for j=1:size(bp,1)        
        if R(j,i)<1
            secc(j,i)=1/sqrt(1-R(j,i));
            out=out+bp(j,i)^2/sqrt(1-R(j,i));
        end
    end
end
%Multiply by dt and double for both hemispheres. For cos square patchbeam, reverse side is 0. 
out=out*(2/size(bp,1)*2/size(bp,2)); %Multiply by dt
%% Using trapezoidal integration
% out=sum(trapz(bp.^2.*secc)*(2/size(bp,1)*2/size(bp,2)));
%% Directivity in dB
out=10*log10(4*pi*m^2/out);
end