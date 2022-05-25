function v=dxyb(u,h)

v = (u([end 1:end-1],[end 1:end-1]) + u(:,:) - u([end 1:end-1],:)-u(:,[end 1:end-1]))/h/h;
