
function v=dyb(u,h)

v = (u(:,:) - u(:,[end 1:end-1]))/h;
