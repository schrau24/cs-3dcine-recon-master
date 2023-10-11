%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A function to do projection normalization and smoothing

%Li Feng, NYU, 12/18/2017

function ZIP=ProjNorm(ZIP);

for ii=1:size(ZIP,3)
    for jj=1:size(ZIP,2)
        tmp=ZIP(:,jj,ii);
        maxprof=max(tmp);
        minprof=min(tmp);
        ZIP(:,jj,ii)=(tmp-minprof)./(maxprof-minprof);
    end
end
for ii=1:size(ZIP,3)
    for jj=1:size(ZIP,1);
        ZIP(jj,:,ii)=smooth(squeeze(ZIP(jj,:,ii)),5,'lowess');
    end
end