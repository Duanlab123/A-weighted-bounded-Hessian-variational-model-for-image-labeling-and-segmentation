
function v=dxb(u,h)

v = (u(:,:) - u([end 1:end-1],:))/h;
