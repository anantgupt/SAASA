function [cov_avg2,cov_min2,cov_max2] = collect_stats(cov_all,range,cfg,Nsnr,Nalg,K_DoA)
cov_avg2=5*log10(squeeze(mean(mean(mean(cov_all(:,range,:,cfg,:,:,:),3),7),2)));
for id=1:K_DoA
    cov_min2(1:Nsnr,1,1:Nalg,id)=5*log10(min(vec(cov_all(:,range,:,cfg,:,id,:))));
    cov_max2(1:Nsnr,1,1:Nalg,id)=5*log10(max(vec(cov_all(:,range,:,cfg,:,id,:))));
end
end