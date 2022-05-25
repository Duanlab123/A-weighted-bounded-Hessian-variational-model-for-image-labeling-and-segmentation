clear all
close all
% profile on
addpath('img'); 
addpath('util');

j = 1;
%%%%%%%%%%%%%%  Data parameters
lambda = 10;
r1 = 1.0;          
r2 = 2.0;
r3 = 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scale = 1;
num_iter = 400;
threshold_res = 1e-3;

Im0 = imread('sh-med-bandages-082404-01.png');
[Ny,Nx,Nc] = size(Im0); 
if Nc>1; Im0=rgb2gray(Im0); end; % Convert color image to gray-scale image
Im0 = double(Im0);

u_ori = Im0/255;

% th=ThdKmeans(u_ori,2);
% th
% 
% c1 = mean(mean(u_ori(u_ori>=th)))
% c2 = mean(mean(u_ori(u_ori<th)))

[l1,l2] = size(u_ori);
Area = l1*l2;

%%%%%%%%%%%%%%  Initialization

u0 = u_ori;

u = rand(l1,l2);
v = u;

scale1 = 5;
ux = dxf(Im0,scale1);
uy = dyf(Im0,scale1);
bta = 1./ sqrt(1 + ux.^2 + uy.^2);

afax = dxf(bta,scale1);
afay = dyf(bta,scale1);
afa = sqrt(afax.^2 + afay.^2);


lambda11 = u0*0; lambda12 = u0*0;
lambda21 = u0*0; lambda22 = u0*0; lambda23 = u0*0; 
lambda3 = u0*0;

x1 = u0*0; x2 = u0*0; 
y1 = u0*0; y2 = u0*0; y3 = u0*0;

scale2 = scale^2;
DtD = abs(psf2otf([1,-1],[l1, l2])).^2 + abs(psf2otf([1;-1],[l1, l2])).^2;
A = (r1/scale2)*DtD + (r2/scale2)*(DtD.^2) + r3;

ave_c = zeros(2,num_iter);

residual = zeros(3,num_iter);
rel_error_L = zeros(3,num_iter);
rel_error_u = zeros(1,num_iter);
energy = zeros(1,num_iter);
rel_error_L(1:3,1) = 1;
rel_error_u(1) = 1;

real_iter = num_iter;
record_flag = 0;

lambda1_L1 = 1;
lambda2_L1 = 1;
lambda3_L1 = 1;
lambda4_L1 = 1;
u_L1 = 1;
% t = cputime;

c1 = 0.35;
c2 = 0.55;
% % % % % % c2 = 0.1788;
% c2 =  0.1106;
tic;
for iter=1:num_iter
    
    %%%%%%%%%%%%%%%%%%%  For c1 and c2
%     [c1,c2] = Cal_Averages(u_ori,v);
%     ave_c(1,iter) = c1;
%     ave_c(2,iter) = c2;    

    %%%%%%%%%%%%%%%%%%%%%  For u
    u_old = u;
    
    g = - dxb(r1*x1+lambda11,scale) -dyb(r1*x2+lambda12,scale) + dxxb(r2*y1+lambda21,scale) + 2*dxyb(r2*y2+lambda22,scale) + dyyb(r2*y3+lambda23,scale) + r3*v +lambda3;   
    g = fftn(g);
    u = real(ifftn(g./A));
       
    %%%%%%%%%%%%%%%%%%%%%  For v
    temp_v = u-( lambda*((u_ori-c1).^2-(u_ori-c2).^2) +lambda3 )/r3;
    v = CV_shrinkage_v(temp_v);  
    
    %%%%%%%%%%%%%%%%%%%%%  For x  
    xx = dxf(u,scale) - lambda11/r1;
    xy = dyf(u,scale) - lambda12/r1;
    x = sqrt(xx.^2 + xy.^2);
    x(x==0) = 1;
    x = max(x - afa/r1,0)./x;
    x1 = xx.*x;
    x2 = xy.*x;
    
    %%%%%%%%%%%%%%%%%%%%%  For y 
    yx = dxxf(u,scale) - lambda21/r2;
    yy = dxyf(u,scale) - lambda22/r2;
    yz = dyyf(u,scale) - lambda23/r2;
    y = sqrt(yx.^2 + 2*yy.^2 + yz.^2);
    y(y==0) = 1;
    y = max(y - bta/r2,0)./y;
    y1 = yx.*y;
    y2 = yy.*y;
    y3 = yz.*y;
    
    %%%%%%%%%%%%%%%%%%%%%  For Lambda
    R11 = x1 - dxf(u,scale);
    R12 = x2 - dyf(u,scale);
    
    R21 = y1 - dxxf(u,scale);
    R22 = y2 - dxyf(u,scale);
    R23 = y3 - dyyf(u,scale);
    
    R3  = v - u;    
    % update
    lambda11_old = lambda11;
    lambda12_old = lambda12;
    lambda21_old = lambda21;
    lambda22_old = lambda22;
    lambda23_old = lambda23;   
    lambda3_old = lambda3; 
    
    lambda11 = lambda11 + r1*R11;
    lambda12 = lambda12 + r1*R12;
   
    lambda21 = lambda21 + r2*R21;
    lambda22 = lambda22 + r2*R22;
    lambda23 = lambda23 + r2*R23;
    
    lambda3  = lambda3  + r3*R3;
    
    %%%%%%%%%%%%%%%%%%%%%  For residual  
    residual(1,iter) = sum( abs(R11(:))+abs(R12(:)) )/Area;
    residual(2,iter) = sum( abs(R21(:))+2*abs(R22(:))+abs(R23(:)) )/Area;
    residual(3,iter) = sum( abs(R3(:)) )/Area;
% 
    %%%%%%%%%%%%%%%%%%%%%  For relative error
    rel_error_u(iter) = sum(sum( abs(u-u_old) ))/u_L1;
    
    rel_error_L(1,iter) = r1*residual(1,iter)/lambda1_L1;
    rel_error_L(2,iter) = r2*residual(2,iter)/lambda2_L1;
    rel_error_L(3,iter) = r3*residual(3,iter)/lambda3_L1;
    
    max_rel_error = max( residual(:,iter) );
    if( max_rel_error < threshold_res )
        if( record_flag==0 )
            real_iter = iter;
            record_flag = 1;
        end
    end
        % L1 norm of Lagrange multipliers.
    if( iter>3 )
        lambda1_L1 = sum( abs(lambda11(:))+abs(lambda12(:)) )/Area;
        lambda2_L1 = sum( abs(lambda21(:))+2*abs(lambda22(:))+ abs(lambda23(:)))/Area;
        lambda3_L1 = sum( abs(lambda3(:)) )/Area;
        u_L1 = sum( abs(u(:)) );
    end
    rel_error_u(iter) = sum(sum( abs(u-u_old) ))/Area;

    energy(iter) = cal_energy(u_ori,c1,c2,v,x1,x2,y1,y2,y3,afa,bta,lambda,scale);
    
%     error(1,iter) = sum(sum(abs(lambda11-lambda11_old)+abs(lambda12-lambda12_old)))/Area;
%     error(2,iter) = sum(sum(abs(lambda21-lambda21_old)+2*abs(lambda22-lambda22_old)+abs(lambda23-lambda23_old)))/Area;
%     error(3,iter) = sum(sum(abs(lambda3-lambda3_old)))/Area;


end
toc
%%%%%%%%%%%%%%%%%% showing results..................
% fprintf(' The iteration number is: %10d  relative error: %15.10f\n ', iter,rel_error_u(iter));
fprintf(' The iteration number is: %10d  relative error: %15.10f\n', real_iter,max(residual(:,iter)) );


xx = 1:iter;
yy = 10:iter;

figure;
plot(yy,log(rel_error_u(10:end) ),'LineWidth',2);
% title('rel error in u');
xlabel('Iteration')
ylabel('Relative error in u^k')
set(gca,'FontWeight','bold')

figure;
plot(yy,log(residual(1,10:end)),'r',yy,log(residual(2,10:end)),'g',yy,log(residual(3,10:end)),'b','LineWidth',2);
% title('rel residuals');
xlabel('Iteration')
ylabel('Relative residuals')
legend('R1','R2','R3','Location','NorthEast')
set(gca,'FontWeight','bold')

figure;
plot(yy,log(rel_error_L(1,10:end)),'r',yy,log(rel_error_L(2,10:end)),'g',yy,log(rel_error_L(3,10:end)),'b','LineWidth',2);
% title('rel errors in L');
xlabel('Iteration')
ylabel('Relative errors in multipliers')
legend('L1','L2','L3','Location','NorthEast')
set(gca,'FontWeight','bold')

figure;
plot(xx,log(energy),'LineWidth',2);
% title('The energy');
xlabel('Iteration')
ylabel('Energy')
set(gca,'FontWeight','bold')

figure,imshow(Im0,[],'border','tight');
hold on;
temp=u<0.5;
contour(temp,'r');


exact1 = imread('sh-med-bandages-082404-01_gt.png');
exact1 = im2bw(exact1);
seg1 = im2bw(u);
figure,imshow(seg1,[],'border','tight')
% % 
[Precision, Recall, JS, SA, Score, k] = pre_rec(seg1,exact1);
     P(j,1)=lambda;
     P(j,2)=r3;
     P(j,3)=Precision;
     P(j,5)=JS;
     P(j,4)=Recall;
     P(j,6)=SA;
     P(j,7)=Score;
     P(j,8)=k;
     j = j+1;
str = ['P.mat'];
save(str,'P');
xlswrite('result.xls',P)
