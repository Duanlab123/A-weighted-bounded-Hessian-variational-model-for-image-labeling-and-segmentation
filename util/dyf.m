
function v=dyf(u,h)

v = (u(:,[2:end 1]) - u(:,:))/h;
