%This function gives magnitude of beampattern of General Planar array
function [patternmra]=bpmra_lite(x,y,SU,SV,sinthetaq,indx,varargin)
%x,y have coordinated of planar array element
rho=1;
if nargin>6
    rho=varargin{1};
end
pattern=zeros(size(SU));
    for i=1:length(x)
        pattern=pattern+exp(1j*(2*pi*rho*(x(i)*(SU)+y(i)*(SV))-sinthetaq(i)));
    end
patternmra=(pattern);
patternmra(indx)=1e-6;
end