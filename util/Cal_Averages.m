
function [c1,c2] = Cal_Averages(u_ori,v)

count1 = sum(sum(v));
c1 = sum(sum( u_ori.*v ));

count2 = sum(sum(1-v));
c2 = sum(sum( u_ori.*(1-v) ));

c1 = c1/(count1+1e-8);
c2 = c2/(count2+1e-8);


