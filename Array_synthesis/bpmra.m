%This function gives magnitude of beampattern of General Planar array
function patternmra=bpmra(x,y,SU,SV,u0,v0,nbits,indx,varargin)
%x,y have coordinated of planar array element
pattern=zeros(size(SU));
SU0=ones(size(SU))*u0;
SV0=ones(size(SV))*v0;
sinthetaq=(x * u0 +y*v0);
rho=1;
if nargin>8
rho=varargin{1};
end
%rho=1.5;%for 30deg
% rho=1.866;%For 60 deg
if nbits<0
    for i=1:length(x)    
        pattern=pattern+exp(1j*2*pi*rho*(x(i)*(SU-SU0)+y(i)*(SV-SV0)));
    end
else
        levels = 2^nbits;
        qlevels = 2.0*pi / levels; % compute quantization levels
        sinthetaq = round(sinthetaq* levels ) .* qlevels; % vector of possible angles
    for i=1:length(x)
        pattern=pattern+exp(1j*(2*pi*rho*(x(i)*(SU)+y(i)*(SV))-sinthetaq(i)));
    end
end
patternmra=(abs(pattern));
patternmra(indx)=1e-6;
end