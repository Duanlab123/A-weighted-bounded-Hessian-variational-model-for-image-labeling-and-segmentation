function v=dxyf(u,h)

v = (u([2:end 1],[2:end 1]) + u(:,:) - u([2:end 1],:)-u(:,[2:end 1]))/h/h;
