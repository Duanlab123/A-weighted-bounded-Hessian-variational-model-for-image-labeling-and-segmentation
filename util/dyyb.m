function v=dyyb(u,h)

v = (u(:,[2:end 1]) + u(:,[end 1:end-1]) - 2*u(:,:))/h/h;
