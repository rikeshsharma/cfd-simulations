%RIKESH SHARMA
%180606
%ME630A                submission on : 20/11/2021
function turbulence2D()
Lx=1;
Ly=1;
Nx=64;
Ny=64;
nxp2=Nx+2;
nyp2=Ny+2;
dx=(Lx)/(Nx-1);
dy=(Ly)/(Ny-1);
x=(0:dx:Lx);
y=(0:dy:Ly);
tmax=0.2;
dt=0.0001;
Nt=tmax/dt;
%since we have to initialize for psi with random number between -1 and 1
%over the whole domain we use a=-1 b=1 and use psi=a+(b-a)*rand(nxp2,nyp2)
%to initialize psi with random number -1 and 1. Here rand(nxp2,nyp2)
%returns a matrix of size (nxp2 X nyp2) with random number between 0 to 1.
a=-1;  b=1;
psi=a+(b-a)*rand(nxp2,nyp2);
psi(Nx+2,:)=psi(2,:);
psi(1,:)=psi(Nx+1,:);
psi(:,Ny+2)=psi(:,2);
psi(:,1)=psi(:,Ny+1);
psi_old=psi;
%Initializing u0, v0, omega0 for t=0.0 using random value of psi over
% the domain
u=zeros(nxp2,nyp2);
u(2:Nx+1,2:Ny+1)=(psi(2:Nx+1,3:Ny+2)-psi(2:Nx+1,1:Ny))/(2*dy);
u(Nx+2,:)=u(2,:);
u(1,:)=u(Nx+1,:);
u(:,Ny+2)=u(:,2);
u(:,1)=u(:,Ny+1);

v=zeros(nxp2,nyp2);
v(2:Nx+1,2:Ny+1)=-(psi(3:Nx+2,2:Ny+1)-psi(1:Nx,2:Ny+1))/(2*dx);
v(Nx+2,:)=v(2,:);
v(1,:)=v(Nx+1,:);
v(:,Ny+2)=v(:,2);
v(:,1)=v(:,Ny+1);

omega=zeros(nxp2,nyp2);
omega(2:Nx+1,2:Ny+1)=(v(3:Nx+2,2:Ny+1)-v(1:Nx,2:Ny+1))/(2*dx) - (u(2:Nx+1,3:Ny+2)-u(2:Nx+1,1:Ny))/(2*dy);
omega(Nx+2,:)=omega(2,:);
omega(1,:)=omega(Nx+1,:);
omega(:,Ny+2)=omega(:,2);
omega(:,1)=omega(:,Ny+1);

% [x,y]=meshgrid(x,y);
% contour(x,y,omega(2:Nx+1,2:Nx+1));
Re=100;
ctr=1;
[x,y]=meshgrid(x,y);
for n=1:Nt
    clf
    %RKW3
    omega_t=omega;
    K1= -u(2:Nx+1,2:Ny+1).*(omega_t(3:Nx+2,2:Ny+1)-omega(1:Nx,2:Ny+1))/(2*dx)...
        -v(2:Nx+1,2:Ny+1).*(omega_t(2:Nx+1,3:Ny+2)-omega(2:Nx+1,1:Ny))/(2*dy)...
        +1/Re*((omega_t(3:Nx+2,2:Ny+1)-2*omega_t(2:Nx+1,2:Ny+1)+omega_t(1:Nx,2:Ny+1))/(dx*dx)...
        +      (omega_t(2:Nx+1,3:Ny+2)-2*omega_t(2:Nx+1,2:Ny+1)+omega_t(2:Nx+1,1:Ny))/(dy*dy));
    
    omega_t(2:Nx+1,2:Ny+1)= omega_t(2:Nx+1,2:Ny+1) +dt/3*K1;
    K1=-5/9*K1   -u(2:Nx+1,2:Ny+1).*(omega_t(3:Nx+2,2:Ny+1)-omega(1:Nx,2:Ny+1))/(2*dx)...
        -v(2:Nx+1,2:Ny+1).*(omega_t(2:Nx+1,3:Ny+2)-omega(2:Nx+1,1:Ny))/(2*dy)...
        +1/Re*((omega_t(3:Nx+2,2:Ny+1)-2*omega_t(2:Nx+1,2:Ny+1)+omega_t(1:Nx,2:Ny+1))/(dx*dx)...
        +      (omega_t(2:Nx+1,3:Ny+2)-2*omega_t(2:Nx+1,2:Ny+1)+omega_t(2:Nx+1,1:Ny))/(dy*dy));
    
    omega_t(2:Nx+1,2:Ny+1) = omega_t(2:Nx+1,2:Ny+1) +15/16*dt*K1;
    K1 = -153/128 *K1  -u(2:Nx+1,2:Ny+1).*(omega_t(3:Nx+2,2:Ny+1)-omega(1:Nx,2:Ny+1))/(2*dx)...
        -v(2:Nx+1,2:Ny+1).*(omega_t(2:Nx+1,3:Ny+2)-omega(2:Nx+1,1:Ny))/(2*dy)...
        +1/Re*((omega_t(3:Nx+2,2:Ny+1)-2*omega_t(2:Nx+1,2:Ny+1)+omega_t(1:Nx,2:Ny+1))/(dx*dx)...
        +      (omega_t(2:Nx+1,3:Ny+2)-2*omega_t(2:Nx+1,2:Ny+1)+omega_t(2:Nx+1,1:Ny))/(dy*dy));
    
    omega(2:Nx+1,2:Ny+1)=omega_t(2:Nx+1,2:Ny+1)+ 8/15*dt*K1;
    omega(Nx+2,:)=omega(2,:);
    omega(1,:)=omega(Nx+1,:);
    omega(:,Ny+2)=omega(:,2);
    omega(:,1)=omega(:,Ny+1);
    %Gauss Siedel
    error=0.5;
    error_old=0;
    changeError=error-error_old;
    epsilon=1e-5;
    while (changeError>epsilon && ctr<10000)
        psi(Nx+2,:)=psi(2,:);
        psi(1,:)=psi(Nx+1,:);
        psi(:,Ny+2)=psi(:,2);
        psi(:,1)=psi(:,Ny+1);
        for i=2:Nx+1
            for j=2:Ny+1
                psi(i,j)=0.25*(psi_old(i+1,j)+psi_old(i,j+1)...
                                  +psi(i-1,j)+psi(i,j-1)+omega(i,j)*dx*dx);
        
            end
        end
%         psi(2:Nx+1,2:Ny+1)=0.25*(psi_old(3:Nx+2,2:Ny+1)+psi_old(2:Nx+1,3:Ny+2)...
%                                   +psi(1:Nx,2:Ny+1)+psi(2:Nx+1,1:Ny)+omega(2:Nx+1,2:Ny+1)*dx*dx);
%         psi(Nx+2,:)=psi(2,:);
%         psi(1,:)=psi(Nx+1,:);
%         psi(:,Ny+2)=psi(:,2);
%         psi(:,1)=psi(:,Ny+1);
        
        error_old=error;
        errorMatrix=(psi(2:Nx+1,2:Ny+1)-psi_old(2:Nx+1,2:Ny+1));
        error = norm(errorMatrix,2);
        changeError=error-error_old;
        ctr=ctr+1;
        psi_old=psi;                 
    end
    %update BC
    omega(Nx+2,:)=omega(2,:);
    omega(1,:)=omega(Nx+1,:);
    omega(:,Ny+2)=omega(:,2);
    omega(:,1)=omega(:,Ny+1);
    psi(Nx+2,:)=psi(2,:);
    psi(1,:)=psi(Nx+1,:);
    psi(:,Ny+2)=psi(:,2);
    psi(:,1)=psi(:,Ny+1);
    
    %calculate u, v
    u(2:Nx+1,2:Ny+1)=(psi(2:Nx+1,3:Ny+2)-psi(2:Nx+1,1:Ny))/(2*dy);
    u(Nx+2,:)=u(2,:);
    u(1,:)=u(Nx+1,:);
    u(:,Ny+2)=u(:,2);
    u(:,1)=u(:,Ny+1);
    v(2:Nx+1,2:Ny+1)=-(psi(3:Nx+2,2:Ny+1)-psi(1:Nx,2:Ny+1))/(2*dx);
    v(Nx+2,:)=v(2,:);
    v(1,:)=v(Nx+1,:);
    v(:,Ny+2)=v(:,2);
    v(:,1)=v(:,Ny+1);
   % Save images Repeat
    contour(x,y,u(2:Nx+1,2:Nx+1),20);
     xlabel('x');
     ylabel('y');
     title('RIKESH SHARMA 180606','\omega 2D Turbulence for Re=100');
     movieFrame(n)=getframe;
    
end
myWriter=VideoWriter('turbulence2D');
myWriter.FrameRate = 30;
open(myWriter);
writeVideo(myWriter, movieFrame);
close(myWriter);

end




