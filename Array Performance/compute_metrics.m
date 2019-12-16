function [cov_all,covc,covcf2,covcf2var] = compute_metrics(covc,covcf2,covcf2var,u_val2,K_DoA,u0,v0)
L_DoA=size(u_val2,2);
err1=kron(u_val2,ones(1,K_DoA))-kron(ones(1,L_DoA),[u0;v0]);
err2=reshape(sum(err1.^2,1),[K_DoA,L_DoA]);
[~,ind_ord]=min(err2,[],2);%Index of u_val's
eu=(u0-u_val2(1,ind_ord));ev=(v0-u_val2(2,ind_ord));
cov_temp=[eu;ev]*[eu.',ev.'];
% cov_all=diag(cov_temp);
covc=covc+cov_temp/K_DoA;% This is storing error in U & V direction of all K targets
for jd=1:K_DoA
    cov_temp_k=[eu(jd);ev(jd)]*[eu(jd),ev(jd)];
    covcf2(:,:,jd)=covcf2(:,:,jd)+cov_temp_k;
    if 0
        covcf2var(:,:,jd,si,1)=covcf2var(:,:,jd,si,1)+(cov_temp_k).^2;
    else
        covcf2var(1,1,jd)=min(covcf2var(1,1,jd),trace(cov_temp_k));
        covcf2var(2,2,jd)=max(covcf2var(2,2,jd),trace(cov_temp_k));
    end
    cov_all(jd)=trace(cov_temp_k);
end
end