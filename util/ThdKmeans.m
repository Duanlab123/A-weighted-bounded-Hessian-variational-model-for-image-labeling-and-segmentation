function th=ThdKmeans(uu,k)

%the number of parts

%   This program is using the K-means method to compute the (k-1)
%   thresholds automatically.
%    
%   Usage: th=ThdKmeans(uu,k);
%
%   Inputs: 
%       - uu: image
%       - k: the number of phases of uu

%
%   Outputs: 
%       - th: the (k-1) thresholds
% 
%
%   Code by: Xiaohao Cai
%   Last updated: 10/11/2012 






% k-means method to separate the given image into k parts
[xLen,yLen]=size(uu);
idx=kmeans(uu(:),k); 
idx=reshape(idx,[xLen,yLen]);

% compute the mean of each part
mean_u(k)=0;
for i=1:k
    temp=(i==idx).*uu;
    mean_u(i)=sum(temp(:))/(sum(temp(:)>0));
end
mean_u=sort(mean_u);

% compute the (k-1) thresholds
th(k-1)=0;
for i=1:k-1
    th(i)=mean(mean_u(i:i+1));
end