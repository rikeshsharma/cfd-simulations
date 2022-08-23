%Rikesh Sharma
%180606
%ME630A (CFD) | < submission date: 22/09/2021 >
function poisson_eq_solution()

% %For 41*41 GRIDS
% %Initialization of domain, grids, coeffients of equations and Boundary Conditions
% domainX=1;
% domainY=1;
% Nx=41;
% Ny=41;
% nxp2=Nx+2;
% nyp2=Ny+2;
% deltaX=domainX/Nx;
% deltaY=domainY/Ny;
% pc = -2*(1/deltaX^2 +1/deltaY^2);
% px = 1/deltaX^2;
% py = 1/deltaY^2;
% x=(-deltaX/2 : deltaX : domainX+deltaX/2);
% y=(-deltaY/2 : deltaY : domainY+deltaY/2);
% bcx=zeros(2,nyp2); bcy=zeros(nxp2,2);
% bcx(1,:)=0.25.*sinh(-5)+(y-0.5).^2.*sinh(10.*(y-0.5))+1;
% bcx(2,:)=0.25.*sinh(5)+(y-0.5).^2.*sinh(10.*(y-0.5))+exp(2.*y);
% bcy(:,1)=0.25.*sinh(-5)+(x-0.5).^2.*sinh(10.*(x-0.5))+1;
% bcy(:,2)=0.25.*sinh(5)+(x-0.5).^2.*sinh(10.*(x-0.5))+exp(2.*x);
% 
% [x,y]=meshgrid(x,y);
% phi=zeros(nxp2,nyp2); phi_old=zeros(nxp2,nyp2);
% Sphi=2.*sinh(10.*(x-0.5))+40.*(x-0.5).*cosh(10.*(x-0.5))+100.*((x-0.5).^2).*sinh(10.*(x-0.5))+2.*sinh(10.*(y-0.5))+40.*(y-0.5).*cosh(10.*(y-0.5))+100.*((y-0.5).^2).*sinh(10.*(y-0.5))+4.*(x.^2+y.^2).*exp(2.*x.*y);
% 
% %Solution for 41*41 Grids%
% residual_Jacobi=zeros(2,1);
% iter_Jacobi=zeros(2,1);
% maxError = 1e-6;
% error=1;
% iter=0;
% %Jacobi method%
% while(error>maxError && iter<100000)
%     iter=iter+1;
%     
%     phi(1,2:Ny+1)  = 2*bcx(1,2:Ny+1) - phi(2,2:Ny+1);
%     phi(nxp2,2:Ny+1) = 2*bcx(2,2:Ny+1) - phi(Nx+1,2:Ny+1);
%     phi(2:Nx+1,1)  = 2*bcy(2:Nx+1,1) - phi(2:Nx+1,2);
%     phi(2:Nx+1,nyp2) = 2*bcy(2:Nx+1,2) - phi(2:Nx+1,Ny+1);
%     
%     for i=2:Nx+1
%         for j=2:Ny+1
%             phi(i,j)=(Sphi(i,j)-(px*phi_old(i-1,j)+px*phi_old(i+1,j)+py*phi_old(i,j-1)+py*phi_old(i,j+1)))/pc;
%         end
%     end
% %     disp(iter);
%      
%     errorMatrix=(phi(2:Nx+1,2:Ny+1)-phi_old(2:Nx+1,2:Ny+1));
%     error = norm(errorMatrix,2);
%     residual_Jacobi(iter,1)=error;
%     iter_Jacobi(iter,1)=iter;
%     phi_old=phi;
% %     disp(error);
% end
% x=(deltaX/2 : deltaX : domainX-deltaX/2);
% y=(deltaY/2 : deltaY : domainY-deltaY/2);
% [x,y]=meshgrid(x,y);
% phi_analytic=(x-0.5).^2.*sinh(10.*(x-0.5))+(y-0.5).^2.*sinh(10.*(y-0.5))+exp(2.*x.*y);
% phi_Jacobi=phi(2:Nx+1,2:Ny+1);
% 
% 
% phi=zeros(nxp2,nyp2); phi_old=zeros(nxp2,nyp2);
% residual_Gauss=zeros(2,1);
% iter_Gauss=zeros(2,1);
% error=1;
% iter=0;
% %Gauss Seidel method%
% while(error>maxError && iter<100000)
%     iter=iter+1;
%     
%     
%     phi(1,2:Ny+1)  = 2*bcx(1,2:Ny+1) - phi(2,2:Ny+1);
%     phi(nxp2,2:Ny+1) = 2*bcx(2,2:Ny+1) - phi(Nx+1,2:Ny+1);
%     phi(2:Nx+1,1)  = 2*bcy(2:Nx+1,1) - phi(2:Nx+1,2);
%     phi(2:Nx+1,nyp2) = 2*bcy(2:Nx+1,2) - phi(2:Nx+1,Ny+1);
%     
%     for i=2:Nx+1
%         for j=2:Ny+1
%             phi(i,j)=(Sphi(i,j)-(px*phi(i-1,j)+px*phi_old(i+1,j)+py*phi(i,j-1)+py*phi_old(i,j+1)))/pc;
%         end
%     end
% %     disp(iter);
%         
%     errorMatrix=(phi(2:Nx+1,2:Ny+1)-phi_old(2:Nx+1,2:Ny+1));
%     error = norm(errorMatrix,2);
%     residual_Gauss(iter,1)=error;
%     iter_Gauss(iter,1)=iter;
%     phi_old=phi;
% %     disp(error);
% end
% 
% phi_Gauss=phi(2:Nx+1,2:Ny+1);
% 
% phi=zeros(nxp2,nyp2); phi_old=zeros(nxp2,nyp2);
% residual_SOR=zeros(2,1);
% iter_SOR=zeros(2,1);
% lambda=1.2;
% error=1;
% iter=0;
% %Succesive Over Relaxation method%
% while(error>maxError && iter<100000)
%     iter=iter+1;
%   
%     phi(1,2:Ny+1)  = 2*bcx(1,2:Ny+1) - phi(2,2:Ny+1);
%     phi(nxp2,2:Ny+1) = 2*bcx(2,2:Ny+1) - phi(Nx+1,2:Ny+1);
%     phi(2:Nx+1,1)  = 2*bcy(2:Nx+1,1) - phi(2:Nx+1,2);
%     phi(2:Nx+1,nyp2) = 2*bcy(2:Nx+1,2) - phi(2:Nx+1,Ny+1);
%     
%     for i=2:Nx+1
%         for j=2:Ny+1
%             phi(i,j)=(1-lambda)*phi(i,j)+lambda*(Sphi(i,j)-(px*phi(i-1,j)+px*phi_old(i+1,j)+py*phi(i,j-1)+py*phi_old(i,j+1)))/pc;
%         end
%     end
% %     disp(iter);
%      
%     errorMatrix=(phi(2:Nx+1,2:Ny+1)-phi_old(2:Nx+1,2:Ny+1));
%     error = norm(errorMatrix,2);
%     residual_SOR(iter,1)=error;
%     iter_SOR(iter,1)=iter;
%     phi_old=phi;
% %     disp(error);
% end
% 
% phi_SOR=phi(2:Nx+1,2:Ny+1);
% 
% phi=zeros(nxp2,nyp2); phi_old=zeros(nxp2,nyp2);
% residual_RBPGS=zeros(2,1);
% iter_RBPGS=zeros(2,1);
% error=1;
% iter=0;
% %RBPGS method%
% while(error>maxError && iter<100000)
%     iter=iter+1;
%     
%  
%     phi(1,2:Ny+1)  = 2*bcx(1,2:Ny+1) - phi(2,2:Ny+1);
%     phi(nxp2,2:Ny+1) = 2*bcx(2,2:Ny+1) - phi(Nx+1,2:Ny+1);
%     phi(2:Nx+1,1)  = 2*bcy(2:Nx+1,1) - phi(2:Nx+1,2);
%     phi(2:Nx+1,nyp2) = 2*bcy(2:Nx+1,2) - phi(2:Nx+1,Ny+1);
%     
%     for i=2:Nx+1
%         m = 2+mod(i,2);
%         for j=m:2:Ny+1
%             phi(i,j) = (Sphi(i,j)-(px*phi_old(i-1,j)+px*phi_old(i+1,j)+py*phi_old(i,j-1)+py*phi_old(i,j+1)))/pc;
%         end
%     end
%     
%     for i=2:Nx+1
%         m = 2+mod(i+1,2);
%         for j=m:2:Ny+1
%             phi(i,j) = (Sphi(i,j)-(px*phi(i-1,j)+px*phi(i+1,j)+py*phi(i,j-1)+py*phi(i,j+1)))/pc;
%         end
%     end
%     
% %     disp(iter);
%      
%     errorMatrix=(phi(2:Nx+1,2:Ny+1)-phi_old(2:Nx+1,2:Ny+1));
%     error = norm(errorMatrix,2);
%     residual_RBPGS(iter,1)=error;
%     iter_RBPGS(iter,1)=iter;
%     phi_old=phi;
% %     disp(error);
% end
% 
% phi_RBPGS=phi(2:Nx+1,2:Ny+1);
% 
% %ADI method%
% phi=zeros(nxp2,nyp2); phi_old=zeros(nxp2,nyp2);
% residual_ADI=zeros(2,1);
% iter_ADI=zeros(2,1);
% error=1;
% iter=0;
% phi_old(1,2:Ny+1)  = 2*bcx(1,2:Ny+1) - phi(2,2:Ny+1);
% phi_old(nxp2,2:Ny+1) = 2*bcx(2,2:Ny+1) - phi(Nx+1,2:Ny+1);
% phi_old(2:Nx+1,1)  = 2*bcy(2:Nx+1,1) - phi(2:Nx+1,2);
% phi_old(2:Nx+1,nyp2) = 2*bcy(2:Nx+1,2) - phi(2:Nx+1,Ny+1);
% 
% ar=zeros(Nx,1);  ar(2:Nx,1)=-px;
% br=zeros(Nx,1);  br(:,1)=-pc;
% cr=zeros(Nx,1);  cr(2:Nx,1)=-px;
% dr=zeros(Nx,1);
% ac=zeros(Ny,1);  ac(2:Ny,1)=-py;
% bc=zeros(Ny,1);  bc(:,1)=-pc;
% cc=zeros(Ny,1);  cc(2:Ny,1)=-py;
% dc=zeros(Ny,1);
% 
% while(error>maxError && iter<100000)
%     iter=iter+1;
% 
%     phi(1,2:Ny+1)  = 2*bcx(1,2:Ny+1) - phi(2,2:Ny+1);
%     phi(nxp2,2:Ny+1) = 2*bcx(2,2:Ny+1) - phi(Nx+1,2:Ny+1);
%     phi(2:Nx+1,1)  = 2*bcy(2:Nx+1,1) - phi(2:Nx+1,2);
%     phi(2:Nx+1,nyp2) = 2*bcy(2:Nx+1,2) - phi(2:Nx+1,Ny+1);
% 
%     for j=2:Ny+1
%         dr(1,1)=py.*phi(2,j-1)+py.*phi_old(2,j+1)-Sphi(2,j)+px*phi(1,j);
%         dr(2:Nx-1,1)=py.*phi(3:Nx,j-1)+py.*phi_old(3:Nx,j+1)-Sphi(3:Nx,j);
%         dr(Nx,1)=py.*phi(Nx+1,j-1)+py.*phi_old(Nx+1,j+1)-Sphi(Nx+1,j)+px*phi(Nx+2,j);
%         phi(2:Nx+1,j)=TDMA(ar,br,cr,dr);
%     end
%     
%     
%     phi(1,2:Ny+1)  = 2*bcx(1,2:Ny+1) - phi(2,2:Ny+1);
%     phi(nxp2,2:Ny+1) = 2*bcx(2,2:Ny+1) - phi(Nx+1,2:Ny+1);
%     phi(2:Nx+1,1)  = 2*bcy(2:Nx+1,1) - phi(2:Nx+1,2);
%     phi(2:Nx+1,nyp2) = 2*bcy(2:Nx+1,2) - phi(2:Nx+1,Ny+1);
%     
%     for i=2:Nx+1
%         dc(1,1)=px.*phi(i-1,2)+px.*phi(i+1,2)-Sphi(i,2)+py*phi(i,1);
%         dc(2:Ny-1,1)=px.*phi(i-1,3:Ny)+px.*phi(i+1,3:Ny)-Sphi(i,3:Ny);
%         dc(Ny,1)=px.*phi(i-1,Ny+1)+px.*phi(i+1,Ny+1)-Sphi(i,Ny+1) + py*phi(i,Ny+2);
%         phi(i,2:Ny+1)=TDMA(ac,bc,cc,dc);
%     end 
%     phi(1,2:Ny+1)  = 2*bcx(1,2:Ny+1) - phi(2,2:Ny+1);
%     phi(nxp2,2:Ny+1) = 2*bcx(2,2:Ny+1) - phi(Nx+1,2:Ny+1);
%     phi(2:Nx+1,1)  = 2*bcy(2:Nx+1,1) - phi(2:Nx+1,2);
%     phi(2:Nx+1,nyp2) = 2*bcy(2:Nx+1,2) - phi(2:Nx+1,Ny+1);
%     
% %     disp(iter);
%      
%     errorMatrix=(phi(2:Nx+1,2:Ny+1)-phi_old(2:Nx+1,2:Ny+1));
%     error = norm(errorMatrix,2);
%     residual_ADI(iter,1)=error;
%     iter_ADI(iter,1)=iter;
%     phi_old=phi;
% %     disp(error);
% end
% 
% phi_ADI=phi(2:Nx+1,2:Ny+1);
% 
% %Plots for 41*41 Grids%
% x=(deltaX/2 : deltaX : domainX-deltaX/2);
% y=(deltaY/2 : deltaY : domainY-deltaY/2);
% [x,y]=meshgrid(x,y);
% 
% %Plotting phi contours
% figure
% subplot(3,2,1)
% contour(x,y,phi_analytic,25);
% xlabel('x');
% ylabel('y');
% title('{\phi}_{analytical}');
% subplot(3,2,2)
% contour(x,y,phi_Jacobi,25);
% xlabel('x');
% ylabel('y');
% title('{\phi} Jacobi (41*41)grids');
% subplot(3,2,3)
% contour(x,y,phi_Gauss,25);
% xlabel('x');
% ylabel('y');
% title('{\phi} Gauss Seidel 41*41)grids');
% subplot(3,2,4)
% contour(x,y,phi_SOR,25);
% xlabel('x');
% ylabel('y');
% title('{\phi} SOR for \lambda = 1.2  (41*41)grids');
% subplot(3,2,5)
% contour(x,y,phi_RBPGS,25);
% xlabel('x');
% ylabel('y');
% title('{\phi} for RBPGS (41*41)grids');
% subplot(3,2,6)
% contour(x,y,phi_ADI,25);
% xlabel('x');
% ylabel('y');
% title('{\phi} for ADI (41*41)grids');
% sgtitle('RIKESH SHARMA 180606');
% 
% %Plotting phi surfaces
% figure
% subplot(3,2,1)
% surfc(x,y,phi_analytic);
% xlabel('x');
% ylabel('y');
% title('{\phi}_{analytical}');
% subplot(3,2,2)
% surfc(x,y,phi_Jacobi);
% xlabel('x');
% ylabel('y');
% title('{\phi} Jacobi (41*41)grids');
% subplot(3,2,3)
% surfc(x,y,phi_Gauss);
% xlabel('x');
% ylabel('y');
% title('{\phi} Gauss Seidel 41*41)grids');
% subplot(3,2,4)
% surfc(x,y,phi_SOR);
% xlabel('x');
% ylabel('y');
% title('{\phi} SOR for \lambda = 1.2  (41*41)grids');
% subplot(3,2,5)
% surfc(x,y,phi_RBPGS);
% xlabel('x');
% ylabel('y');
% title('{\phi} for RBPGS (41*41)grids');
% subplot(3,2,6)
% surfc(x,y,phi_ADI);
% xlabel('x');
% ylabel('y');
% title('{\phi} for ADI (41*41)grids');
% sgtitle('RIKESH SHARMA 180606');
% 
% 
% %Plotting Errors
% figure
% subplot(3,2,1)
% contour(x,y,phi_Jacobi-phi_analytic,'ShowText','on');
% xlabel('x');
% ylabel('y');
% title('Error for Jacobi (41*41)grids');
% subplot(3,2,2)
% contour(x,y,phi_Gauss-phi_analytic,'ShowText','on');
% xlabel('x');
% ylabel('y');
% title('Error for Gauss Seidel (41*41)grids');
% subplot(3,2,3)
% contour(x,y,phi_SOR-phi_analytic,'ShowText','on');
% xlabel('x');
% ylabel('y');
% title('Error for SOR (\lambda=1.2) (41*41)grids');
% subplot(3,2,4)
% contour(x,y,phi_RBPGS-phi_analytic,'ShowText','on');
% xlabel('x');
% ylabel('y');
% title('Error for RBPGS (41*41)grids');
% subplot(3,2,5)
% x=(9*deltaX/2 : deltaX : domainX-deltaX/2);
% y=(9*deltaY/2 : deltaY : domainY-deltaY/2);
% [x,y]=meshgrid(x,y);
% contour(x,y,phi_ADI(5:Nx,5:Ny)-phi_analytic(5:Nx,5:Ny),'ShowText','on');
% xlabel('x');
% ylabel('y');
% title('Error for ADI (41*41)grids');
% sgtitle('RIKESH SHARMA 180606');
% 
% 
% %Ploting residuals for rate of convergence
% figure
% semilogy(iter_Jacobi ,residual_Jacobi);
% hold on
% semilogy(iter_Gauss  ,residual_Gauss);
% semilogy(iter_SOR    ,residual_SOR);
% semilogy(iter_RBPGS  ,residual_RBPGS);
% semilogy(iter_ADI    ,residual_ADI);
% yline(1e-6,'-.k','Tolerance = 10^-^6');
% legend('Jacobi','Gauss-Seidel','SOR for \lambda=1.2','RBPGS','ADI','Tolerance line');
% xlabel('Number of iterations');
% ylabel('Residual on Log-scale');
% title('Residual vs Number of iterations for 41*41 grids on semilog plot(Rikesh Sharma 180606)');
% 
% save('phi41.mat','phi_analytic','phi_Jacobi','phi_Gauss','phi_SOR','phi_RBPGS','phi_ADI');
% save('residual41.mat','residual_Jacobi','residual_Gauss','residual_SOR','residual_RBPGS','residual_ADI');
% %Solution for 41*41 Grids Ends%

%For 81*81 GRIDS
%Initialization of domain, grids, coeffients of equations and Boundary Conditions
domainX=1;
domainY=1;
Nx=200;
Ny=200;
nxp2=Nx+2;
nyp2=Ny+2;
deltaX=domainX/Nx;
deltaY=domainY/Ny;
pc = -2*(1/deltaX^2 +1/deltaY^2);
px = 1/deltaX^2;
py = 1/deltaY^2;
x=(-deltaX/2 : deltaX : domainX+deltaX/2);
y=(-deltaY/2 : deltaY : domainY+deltaY/2);
bcx=zeros(2,nyp2); bcy=zeros(nxp2,2);
bcx(1,:)=0.25.*sinh(-5)+(y-0.5).^2.*sinh(10.*(y-0.5))+1;
bcx(2,:)=0.25.*sinh(5)+(y-0.5).^2.*sinh(10.*(y-0.5))+exp(2.*y);
bcy(:,1)=0.25.*sinh(-5)+(x-0.5).^2.*sinh(10.*(x-0.5))+1;
bcy(:,2)=0.25.*sinh(5)+(x-0.5).^2.*sinh(10.*(x-0.5))+exp(2.*x);

[x,y]=meshgrid(x,y);
phi=zeros(nxp2,nyp2); phi_old=zeros(nxp2,nyp2);
Sphi=2.*sinh(10.*(x-0.5))+40.*(x-0.5).*cosh(10.*(x-0.5))+100.*((x-0.5).^2).*sinh(10.*(x-0.5))+2.*sinh(10.*(y-0.5))+40.*(y-0.5).*cosh(10.*(y-0.5))+100.*((y-0.5).^2).*sinh(10.*(y-0.5))+4.*(x.^2+y.^2).*exp(2.*x.*y);

%Solution for 41*41 Grids%
residual_Jacobi=zeros(2,1);
iter_Jacobi=zeros(2,1);
maxError = 1e-6;
error=1;
iter=0;
%Jacobi method%
while(error>maxError && iter<100000)
    iter=iter+1;
    
    phi(1,2:Ny+1)  = 2*bcx(1,2:Ny+1) - phi(2,2:Ny+1);
    phi(nxp2,2:Ny+1) = 2*bcx(2,2:Ny+1) - phi(Nx+1,2:Ny+1);
    phi(2:Nx+1,1)  = 2*bcy(2:Nx+1,1) - phi(2:Nx+1,2);
    phi(2:Nx+1,nyp2) = 2*bcy(2:Nx+1,2) - phi(2:Nx+1,Ny+1);
    
    for i=2:Nx+1
        for j=2:Ny+1
            phi(i,j)=(Sphi(i,j)-(px*phi_old(i-1,j)+px*phi_old(i+1,j)+py*phi_old(i,j-1)+py*phi_old(i,j+1)))/pc;
        end
    end
%     disp(iter);
     
    errorMatrix=(phi(2:Nx+1,2:Ny+1)-phi_old(2:Nx+1,2:Ny+1));
    error = norm(errorMatrix,2);
    residual_Jacobi(iter,1)=error;
    iter_Jacobi(iter,1)=iter;
    phi_old=phi;
%     disp(error);
end
x=(deltaX/2 : deltaX : domainX-deltaX/2);
y=(deltaY/2 : deltaY : domainY-deltaY/2);
[x,y]=meshgrid(x,y);
phi_analytic=(x-0.5).^2.*sinh(10.*(x-0.5))+(y-0.5).^2.*sinh(10.*(y-0.5))+exp(2.*x.*y);
phi_Jacobi=phi(2:Nx+1,2:Ny+1);

phi=zeros(nxp2,nyp2); phi_old=zeros(nxp2,nyp2);
residual_Gauss=zeros(2,1);
iter_Gauss=zeros(2,1);
error=1;
iter=0;
%Gauss Seidel method%
while(error>maxError && iter<100000)
    iter=iter+1;
    
    
    phi(1,2:Ny+1)  = 2*bcx(1,2:Ny+1) - phi(2,2:Ny+1);
    phi(nxp2,2:Ny+1) = 2*bcx(2,2:Ny+1) - phi(Nx+1,2:Ny+1);
    phi(2:Nx+1,1)  = 2*bcy(2:Nx+1,1) - phi(2:Nx+1,2);
    phi(2:Nx+1,nyp2) = 2*bcy(2:Nx+1,2) - phi(2:Nx+1,Ny+1);
    
    for i=2:Nx+1
        for j=2:Ny+1
            phi(i,j)=(Sphi(i,j)-(px*phi(i-1,j)+px*phi_old(i+1,j)+py*phi(i,j-1)+py*phi_old(i,j+1)))/pc;
        end
    end
%     disp(iter);
        
    errorMatrix=(phi(2:Nx+1,2:Ny+1)-phi_old(2:Nx+1,2:Ny+1));
    error = norm(errorMatrix,2);
    residual_Gauss(iter,1)=error;
    iter_Gauss(iter,1)=iter;
    phi_old=phi;
%     disp(error);
end

phi_Gauss=phi(2:Nx+1,2:Ny+1);

phi=zeros(nxp2,nyp2); phi_old=zeros(nxp2,nyp2);
residual_SOR=zeros(2,1);
iter_SOR=zeros(2,1);
lambda=1.2;
error=1;
iter=0;
%Succesive Over Relaxation method%
while(error>maxError && iter<100000)
    iter=iter+1;
  
    phi(1,2:Ny+1)  = 2*bcx(1,2:Ny+1) - phi(2,2:Ny+1);
    phi(nxp2,2:Ny+1) = 2*bcx(2,2:Ny+1) - phi(Nx+1,2:Ny+1);
    phi(2:Nx+1,1)  = 2*bcy(2:Nx+1,1) - phi(2:Nx+1,2);
    phi(2:Nx+1,nyp2) = 2*bcy(2:Nx+1,2) - phi(2:Nx+1,Ny+1);
    
    for i=2:Nx+1
        for j=2:Ny+1
            phi(i,j)=(1-lambda)*phi(i,j)+lambda*(Sphi(i,j)-(px*phi(i-1,j)+px*phi_old(i+1,j)+py*phi(i,j-1)+py*phi_old(i,j+1)))/pc;
        end
    end
%     disp(iter);
     
    errorMatrix=(phi(2:Nx+1,2:Ny+1)-phi_old(2:Nx+1,2:Ny+1));
    error = norm(errorMatrix,2);
    residual_SOR(iter,1)=error;
    iter_SOR(iter,1)=iter;
    phi_old=phi;
%     disp(error);
end

phi_SOR=phi(2:Nx+1,2:Ny+1);

phi=zeros(nxp2,nyp2); phi_old=zeros(nxp2,nyp2);
residual_RBPGS=zeros(2,1);
iter_RBPGS=zeros(2,1);
error=1;
iter=0;
%RBPGS method%
while(error>maxError && iter<100000)
    iter=iter+1;
    
 
    phi(1,2:Ny+1)  = 2*bcx(1,2:Ny+1) - phi(2,2:Ny+1);
    phi(nxp2,2:Ny+1) = 2*bcx(2,2:Ny+1) - phi(Nx+1,2:Ny+1);
    phi(2:Nx+1,1)  = 2*bcy(2:Nx+1,1) - phi(2:Nx+1,2);
    phi(2:Nx+1,nyp2) = 2*bcy(2:Nx+1,2) - phi(2:Nx+1,Ny+1);
    
    for i=2:Nx+1
        m = 2+mod(i,2);
        for j=m:2:Ny+1
            phi(i,j) = (Sphi(i,j)-(px*phi_old(i-1,j)+px*phi_old(i+1,j)+py*phi_old(i,j-1)+py*phi_old(i,j+1)))/pc;
        end
    end
    
    for i=2:Nx+1
        m = 2+mod(i+1,2);
        for j=m:2:Ny+1
            phi(i,j) = (Sphi(i,j)-(px*phi(i-1,j)+px*phi(i+1,j)+py*phi(i,j-1)+py*phi(i,j+1)))/pc;
        end
    end
    
%     disp(iter);
     
    errorMatrix=(phi(2:Nx+1,2:Ny+1)-phi_old(2:Nx+1,2:Ny+1));
    error = norm(errorMatrix,2);
    residual_RBPGS(iter,1)=error;
    iter_RBPGS(iter,1)=iter;
    phi_old=phi;
%     disp(error);
end

phi_RBPGS=phi(2:Nx+1,2:Ny+1);

%ADI method%
phi=zeros(nxp2,nyp2); phi_old=zeros(nxp2,nyp2);
residual_ADI=zeros(2,1);
iter_ADI=zeros(2,1);
error=1;
iter=0;
phi_old(1,2:Ny+1)  = 2*bcx(1,2:Ny+1) - phi(2,2:Ny+1);
phi_old(nxp2,2:Ny+1) = 2*bcx(2,2:Ny+1) - phi(Nx+1,2:Ny+1);
phi_old(2:Nx+1,1)  = 2*bcy(2:Nx+1,1) - phi(2:Nx+1,2);
phi_old(2:Nx+1,nyp2) = 2*bcy(2:Nx+1,2) - phi(2:Nx+1,Ny+1);

ar=zeros(Nx,1);  ar(2:Nx,1)=-px;
br=zeros(Nx,1);  br(:,1)=-pc;
cr=zeros(Nx,1);  cr(2:Nx,1)=-px;
dr=zeros(Nx,1);
ac=zeros(Ny,1);  ac(2:Ny,1)=-py;
bc=zeros(Ny,1);  bc(:,1)=-pc;
cc=zeros(Ny,1);  cc(2:Ny,1)=-py;
dc=zeros(Ny,1);

while(error>maxError && iter<100000)
    iter=iter+1;

    phi(1,2:Ny+1)  = 2*bcx(1,2:Ny+1) - phi(2,2:Ny+1);
    phi(nxp2,2:Ny+1) = 2*bcx(2,2:Ny+1) - phi(Nx+1,2:Ny+1);
    phi(2:Nx+1,1)  = 2*bcy(2:Nx+1,1) - phi(2:Nx+1,2);
    phi(2:Nx+1,nyp2) = 2*bcy(2:Nx+1,2) - phi(2:Nx+1,Ny+1);

    for j=2:Ny+1
        dr(1,1)=py.*phi(2,j-1)+py.*phi_old(2,j+1)-Sphi(2,j)+px*phi(1,j);
        dr(2:Nx-1,1)=py.*phi(3:Nx,j-1)+py.*phi_old(3:Nx,j+1)-Sphi(3:Nx,j);
        dr(Nx,1)=py.*phi(Nx+1,j-1)+py.*phi_old(Nx+1,j+1)-Sphi(Nx+1,j)+px*phi(Nx+2,j);
        phi(2:Nx+1,j)=TDMA(ar,br,cr,dr);
    end
    
    
    phi(1,2:Ny+1)  = 2*bcx(1,2:Ny+1) - phi(2,2:Ny+1);
    phi(nxp2,2:Ny+1) = 2*bcx(2,2:Ny+1) - phi(Nx+1,2:Ny+1);
    phi(2:Nx+1,1)  = 2*bcy(2:Nx+1,1) - phi(2:Nx+1,2);
    phi(2:Nx+1,nyp2) = 2*bcy(2:Nx+1,2) - phi(2:Nx+1,Ny+1);
    
    for i=2:Nx+1
        dc(1,1)=px.*phi(i-1,2)+px.*phi(i+1,2)-Sphi(i,2)+py*phi(i,1);
        dc(2:Ny-1,1)=px.*phi(i-1,3:Ny)+px.*phi(i+1,3:Ny)-Sphi(i,3:Ny);
        dc(Ny,1)=px.*phi(i-1,Ny+1)+px.*phi(i+1,Ny+1)-Sphi(i,Ny+1) + py*phi(i,Ny+2);
        phi(i,2:Ny+1)=TDMA(ac,bc,cc,dc);
    end 
    phi(1,2:Ny+1)  = 2*bcx(1,2:Ny+1) - phi(2,2:Ny+1);
    phi(nxp2,2:Ny+1) = 2*bcx(2,2:Ny+1) - phi(Nx+1,2:Ny+1);
    phi(2:Nx+1,1)  = 2*bcy(2:Nx+1,1) - phi(2:Nx+1,2);
    phi(2:Nx+1,nyp2) = 2*bcy(2:Nx+1,2) - phi(2:Nx+1,Ny+1);
    
%     disp(iter);
     
    errorMatrix=(phi(2:Nx+1,2:Ny+1)-phi_old(2:Nx+1,2:Ny+1));
    error = norm(errorMatrix,2);
    residual_ADI(iter,1)=error;
    iter_ADI(iter,1)=iter;
    phi_old=phi;
%     disp(error);
end

phi_ADI=phi(2:Nx+1,2:Ny+1);

x=(deltaX/2 : deltaX : domainX-deltaX/2);
y=(deltaY/2 : deltaY : domainY-deltaY/2);
[x,y]=meshgrid(x,y);

%Plots for 81*81 Grids%
%Ploting phi contours
figure
subplot(3,2,1)
contour(x,y,phi_analytic,25);
xlabel('x');
ylabel('y');
title('{\phi}_{analytical}');
subplot(3,2,2)
contour(x,y,phi_Jacobi,25);
xlabel('x');
ylabel('y');
title('{\phi} for Jacobi (81*81)grids');
subplot(3,2,3)
contour(x,y,phi_Gauss,25);
xlabel('x');
ylabel('y');
title('{\phi} for Gauss Seidel (81*81)grids');
subplot(3,2,4)
contour(x,y,phi_SOR,25);
xlabel('x');
ylabel('y');
title('{\phi} for SOR (81*81)grids');
subplot(3,2,5)
contour(x,y,phi_RBPGS,25);
xlabel('x');
ylabel('y');
title('{\phi} for RBPGS (81*81)grids');
subplot(3,2,6)
contour(x,y,phi_ADI,25);
xlabel('x');
ylabel('y');
title('{\phi} for ADI (81*81)grids');
sgtitle('RIKESH SHARMA 180606');


%Ploting phi surfaces
figure
subplot(3,2,1)
surfc(x,y,phi_analytic);
xlabel('x');
ylabel('y');
title('{\phi}_{analytical}');
subplot(3,2,2)
surfc(x,y,phi_Jacobi);
xlabel('x');
ylabel('y');
title('{\phi} for Jacobi (81*81)grids');
subplot(3,2,3)
surfc(x,y,phi_Gauss);
xlabel('x');
ylabel('y');
title('{\phi} for Gauss Seidel (81*81)grids');
subplot(3,2,4)
surfc(x,y,phi_SOR);
xlabel('x');
ylabel('y');
title('{\phi} for SOR (81*81)grids');
subplot(3,2,5)
surfc(x,y,phi_RBPGS);
xlabel('x');
ylabel('y');
title('{\phi} for RBPGS (81*81)grids');
subplot(3,2,6)
surfc(x,y,phi_ADI);
xlabel('x');
ylabel('y');
title('{\phi} for ADI (81*81)grids');
sgtitle('RIKESH SHARMA 180606');

%Ploting Error contours
figure
subplot(3,2,1)
contour(x,y,phi_Jacobi-phi_analytic,'ShowText','on');
xlabel('x');
ylabel('y');
title('Error for Jacobi (81*81)grids');
subplot(3,2,2)
contour(x,y,phi_Gauss-phi_analytic,'ShowText','on');
xlabel('x');
ylabel('y');
title('Error for Gauss Seidel (81*81)grids');
subplot(3,2,3)
contour(x,y,phi_SOR-phi_analytic,'ShowText','on');
xlabel('x');
ylabel('y');
title('Error for SOR (81*81)grids');
subplot(3,2,4)
contour(x,y,phi_RBPGS-phi_analytic,'ShowText','on');
xlabel('x');
ylabel('y');
title('Error for RBPGS (81*81)grids');
subplot(3,2,5)
x=(13*deltaX/2 : deltaX : domainX-deltaX/2);
y=(13*deltaY/2 : deltaY : domainY-deltaY/2);
[x,y]=meshgrid(x,y);
contour(x,y,phi_ADI(7:Nx,7:Ny)-phi_analytic(7:Nx,7:Ny),'ShowText','on');
xlabel('x');
ylabel('y');
title('Error for ADI (81*81)grids');
sgtitle('RIKESH SHARMA 180606');

%Ploting residuals
figure
semilogy(iter_Jacobi ,residual_Jacobi);
hold on
semilogy(iter_Gauss  ,residual_Gauss);
semilogy(iter_SOR    ,residual_SOR);
semilogy(iter_RBPGS  ,residual_RBPGS);
semilogy(iter_ADI    ,residual_ADI);
yline(1e-6,'-.k','Tolerance = 10^-^6');
legend('Jacobi','Gauss-Seidel','SOR for \lambda=1.2','RBPGS','ADI','Tolerance line');
xlabel('Number of iterations');
ylabel('Residual on Log-scale');
title('Residual vs Number of iterations for (81*81 grids) on semilog plot (Rikesh Sharma 180606)');

save('phi81.mat','phi_analytic','phi_Jacobi','phi_Gauss','phi_SOR','phi_RBPGS','phi_ADI');
save('residual81.mat','residual_Jacobi','residual_Gauss','residual_SOR','residual_RBPGS','residual_ADI');

%Solutions for (81*81) Grids ends%


end

    
    
    
    
    
    
    
    
