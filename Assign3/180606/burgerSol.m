function bergerSol()

Re=10;
xmax=1;
nu=1/Re;
dt=0.001;
dx=0.001;
x=-xmax:dx:xmax;
S=nu*dt/(dx*dx);
N=1/dt+1;
gridP=2*xmax/dx +1;
a=zeros(gridP,N);
b=zeros(gridP,N);
c=zeros(gridP,N);
d=zeros(gridP,N);
u=zeros(gridP,N);
u(1:xmax/dx,1)=1.0;
u(xmax/dx+1:gridP,1)=0.0;
u(1,:)=1.0;
u(gridP,:)=0.0;

for i=2:N
    
   b(1)=1; 
%    b(47)=1;
%     a(47)=0;
    c(1)=0.0;
    d(1)=1.0; 
%     d(47)=0.0;
    
    a(2:gridP-1,i-1)=-0.25*(dt/dx)*u(1:gridP-2,i-1)-0.5*S;
    b(2:gridP-1,i-1)=1+S;
    c(2:gridP-1,i-1)=0.25*(dt/dx)*u(3:gridP,i-1)-0.5*S;
    d(2:gridP-1,i-1)=0.5*S*u(1:gridP-2,i-1)+(1-S)*u(2:gridP-1,i-1)+0.5*S*u(3:gridP,i-1);
    
%   d(2,i-1)=0.5*S*u(1,i-1)+(1-S)*u(2,i-1)+0.5*S*u(3,i-1)+(0.25*(dt/dx)*u(1,i-1)-0.5*S)*u(1,i);
%   d(gridP-1,i-1)=0.5*S*u(gridP-2,i-1)+(1-S)*u(gridP-1,i-1)+0.5*S*u(gridP,i-1)-(0.25*(dt/dx)*u(gridP,i-1)-0.5*S)*u(gridP,i);
  
    u(:,i)=TDMA(c(:,i-1),b(:,i-1),a(:,i-1),d(:,i-1));
   
%   u(1:xmax/dx,1)=1.0;
%   u(xmax/dx+1:gridP,1)=0.0;
%   u(1,:)=1.0;
%   u(gridP,:)=0.0;
end
u
plot(x,u(:,2));




end
