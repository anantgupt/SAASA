function [out]=patchbeam(u,v)
% eps=0.001;
% p=atan(v./u);t=asin(sqrt(u.^2+v.^2));
% W=2e-3;
% L=2e-3;
% k=2*pi*60e9/3e8;%Free space Wavenumber
% Et=@(p,t)sin(k*W*sin(t).*sin(p)/2)./(k*W*sin(t)*sin(p)/2)*cos(k*L/2*sin(t)*cos(p))*cos(p);
% Ep=@(p,t)sin(k*W*sin(t).*sin(p)/2)./(k*W*sin(t)*sin(p)/2)*cos(k*L/2*sin(t)*cos(p))*cos(t)*sin(p);
% fmag=sqrt(Et(p,t).^2+Ep(p,t).^2); %Total Field Magnitude
%  if(u.^2+v.^2>1)
%      out=eps;
%  else
%         if t==0 
%             out=1;
%         elseif p==0
%             out=abs(cos(k*L/2*sin(t)));
%         else
%             out=fmag;
%         end
%  end

 out=real(sqrt(ones(size(u))-u.^2-v.^2));%cosine pattern
%  out=real(cos(t));%cosine pattern
end
