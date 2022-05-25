DtD = abs(psf2otf([1,-1],[256 256])).^2 + abs(psf2otf([1;-1],[256 256])).^2;
Z = DtD.^2;
figure;imshow(DtD,[]);
figure;imshow(Z,[])
l1 = 256;
l2 = 256;
z1p = -1 + exp( sqrt(-1)*2*pi*[1:l1]/l1);
z1n =  1 - exp(-sqrt(-1)*2*pi*[1:l1]/l1);
z2p = -1 + exp( sqrt(-1)*2*pi*[1:l2]/l2);
z2n =  1 - exp(-sqrt(-1)*2*pi*[1:l2]/l2);

A = zeros(l1,l2);
for i=1:l1
    for j=1:l2
        A(i,j) = -2*( cos(2*pi*i/l1)+cos(2*pi*j/l2)-2 )+ 4*(r2+mu)*( cos(2*pi*(i-1)/l1)+cos(2*pi*(j-1)/l2)-2 ).^2  + r3;
    end
end
figure;imshow(A,[]);
c = A.^2;
figure;imshow(c,[]);