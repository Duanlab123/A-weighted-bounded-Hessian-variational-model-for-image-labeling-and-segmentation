
function energy = cal_energy(u_ori,c1,c2,v,x1,x2,y1,y2,y3,afa,bta,lambda,scale)

energy = 0;
norm_x = sqrt(x1.^2 + x2.^2);
normx = afa.*norm_x;
norm_y = sqrt(y1.^2 + 2*(y2.^2) + y3.^2);
normy = bta.*norm_y;
tem = lambda*((u_ori-c1).^2.*v + (u_ori-c2).^2.*(1-v)) + normx + normy;
energy = sum(sum(tem));
% energy = energy*scale*scale;


