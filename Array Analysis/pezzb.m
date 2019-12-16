%% Pe for Non-coherent BHT (using magnitudes of S(u)*S(u+d))
function Pe = pezzb(Rd,s)
% RD: Correlation
% s: SNR of complex noise added
a=sqrt((1-sqrt(1-Rd.^2))*s/2);
b=sqrt((1+sqrt(1-Rd.^2))*s/2);
% a=sqrt(Rd.^2*s);
% b=sqrt(s);
Pe=marcumq(a,b,1)-0.5*exp(-(a.^2+b.^2)/2).*besseli(0,a.*b);
end