%RIKESH SHARMA
%180606
%ME630A

function burgerSolutions()
xmax=1;
dx=0.01;
x=-xmax:dx:xmax;

Re=10;
ue10=explicitBurger(Re);
ui10=implicitBurger(Re);
uExact10=exactBurgerSol(Re,0.47);

Re=50;
ue50=explicitBurger(Re);
ui50=implicitBurger(Re);
uExact50=exactBurgerSol(Re,0.47);

figure
plot(x,ue10(:,4700));
hold on
plot(x,ue50(:,4700));
plot(x,ui10(:,47));
plot(x,ui50(:,47));
plot(x,uExact10);
plot(x,uExact50);
xlabel(' X from (-1 to 1) ');
ylabel('U for x at time t=0.47 sec');
title('RIKESH SHARMA 180606','U for x at time t=0.47 sec');
legend('Explicit Re=10(FTFS)','Explicit Re=50(FTFS)','Implicit Re=10 (CN)','Implicit Re=50 (CN)','Exact Re=10','Exact Re=50');
xlim([-1 1]);
ylim([-0.01 1.01]);

figure
plot(x,ue10(:,4700));
hold on
plot(x,ui10(:,47));
plot(x,uExact10);
xlabel(' X from (-1 to 1) ');
ylabel('U for x at time t=0.47 sec');
title('RIKESH SHARMA 180606','U for x at time t=0.47 sec');
legend('Explicit Re=10(FTFS)','Implicit Re=10 (CN)','Exact Re=10');
xlim([-1 1]);
ylim([-0.01 1.01]);

figure
plot(x,ue50(:,4700));
hold on
plot(x,ui50(:,47));
plot(x,uExact50);
xlabel(' X from (-1 to 1) ');
ylabel('U for x at time t=0.47 sec');
title('RIKESH SHARMA 180606','U for x at time t=0.47 sec');
legend('Explicit Re=50(FTFS)','Implicit Re=50 (CN)','Exact Re=50');
xlim([-1 1]);
ylim([-0.01 1.01]);

figure
plot(x,ue10(:,4700));
hold on
plot(x,ue50(:,4700));
xlabel(' X from (-1 to 1) ');
ylabel('U for x at time t=0.47 sec');
title('RIKESH SHARMA 180606','U for x at time t=0.47 sec');
legend('Explicit Re=10(FTFS)','Explicit Re=50(FTFS)');
xlim([-1 1]);
ylim([-0.01 1.01]);

figure
plot(x,ui10(:,47));
hold on
plot(x,ui50(:,47));
xlabel(' X from (-1 to 1) ');
ylabel('U for x at time t=0.47 sec');
title('RIKESH SHARMA 180606','U for x at time t=0.47 sec');
legend('Implicit Re=10 (CN)','Implicit Re=50 (CN)');
xlim([-1 1]);
ylim([-0.01 1.01]);
end

function u=explicitBurger(Re)
xmax=1;
nu=1/Re;
dt=0.0001;
dx=0.01;

s=nu*dt/(dx*dx);

N=(1/dt)+1;
gridP=2*xmax/dx +1;

u=zeros(gridP,N);
u(1:xmax/dx,1)=1.0;
u(xmax/dx+1:gridP,1)=0.0;
u(1,:)=1.0;
u(gridP,:)=0.0;

for i=1:N-1
       
%      u(2:gridP-1,i+1)=u(2:gridP-1,i)...
%          -(dt/dx)*u(2:gridP-1,i).*(u(2:gridP-1,i)-u(1:gridP-2,i))...
%          +0.5*(dt^2/dx^2)*u(2:gridP-1,i+1).^2.*(u(3:gridP,i+1)-2*u(2:gridP-1,i+1)+u(1:gridP-2,i+1))...
%          +s*(u(3:gridP,i)-2*u(2:gridP-1,i)+u(1:gridP-2,i));
    u(2:gridP-1,i+1)=u(2:gridP-1,i)...
        -(dt/dx)*u(2:gridP-1,i).*(u(3:gridP,i)-u(2:gridP-1,i))...
         +0.5*(dt^2/dx^2)*u(2:gridP-1,i+1).^2.*(u(3:gridP,i+1)-2*u(2:gridP-1,i+1)+u(1:gridP-2,i+1))...
         +s*(u(3:gridP,i)-2*u(2:gridP-1,i)+u(1:gridP-2,i));
         
end

end
function u=implicitBurger(Re)
xmax=1;
nu=1/Re;
dt=0.01;
dx=0.01;

s=nu*dt/(dx*dx);

N=(1/dt)+1;
gridP=2*xmax/dx +1;

u=zeros(gridP,N);
u(1:xmax/dx,1)=1.0;
u(xmax/dx+1:gridP,1)=0.0;
u(1,:)=1.0;
u(gridP,:)=0.0;
a=zeros(gridP-2);
b=zeros(gridP-2);
c=zeros(gridP-2);
d=zeros(gridP-2);

for i=1:N-1
    b(1:gridP-2) = ( 1-u(2:gridP-1,i)*dt/(2*dx)+s);
    a(2:gridP-2) = -0.5*s;
    c(1:gridP-3) = (0.5*u(3:gridP-1,i)*(dt/dx)-0.5*s);
    d(1)         = 0.5*s*u(1,i)+(1-s)*u(2,i)+0.5*s*u(3,i)+u(1,i+1)*0.5*s;
    d(gridP-2)   = 0.5*s*u(gridP-2,i)+(1-s)*u(gridP-1,i)+0.5*s*u(gridP,i)- u(gridP,i+1)*(u(gridP,i)*0.5*(dt/dx)-0.5*s);
    d(2:gridP-3) = 0.5*s*u(2:gridP-3,i)+(1-s)*u(3:gridP-2,i)+0.5*s*u(4:gridP-1,i);
    
    u(2:gridP-1,i+1) = TDMA(a,b,c,d);
    
end


end
function uExact=exactBurgerSol(Re,t)

xmax=1;
dx=0.01;
x=-xmax:dx:xmax;
N=2*xmax/dx+1;
uExact=zeros(N,1);
for i = 1:N
    f = @(z)  ((x(i)-z)/t).*exp(-0.5*Re*(z+0.5*((x(i)-z).^2)/t)) ;
    g = @(z)  ((x(i)-z)/t).*exp(-0.5*Re*(0.5*((x(i)-z).^2)/t)) ;
    h = @(z)  exp(-0.5*Re*(z+0.5*((x(i)-z).^2)/t));
    l = @(z)  exp(-0.5*Re*(0.5*((x(i)-z).^2)/t));
    
    fi = integral(f,-inf,0);
    gi = integral(g,0,inf);
    hi = integral(h,-inf,0);
    li = integral(l,0,inf);
    
    uExact(i)=(fi+gi)/(hi+li);
end

end
function x = TDMA(a,b,c,d)
%a, b, c are the column vectors for the compressed tridiagonal matrix, d is the right vector
n = length(b); % n is the number of rows
 
% Modify the first-row coefficients
c(1) = c(1) / b(1);    % Division by zero risk.
d(1) = d(1) / b(1);    % Division by zero would imply a singular matrix.
 
for i = 2:n-1
    mult = b(i) - a(i) * c(i-1);
    c(i) = c(i) / mult;
    d(i) = (d(i) - a(i) * d(i-1)) / mult;
end
 
d(n) = (d(n) - a(n) * d(n-1))/( b(n) - a(n) * c(n-1));
 
% Now back substitute.
x(n) = d(n);
for i = n-1:-1:1
    x(i) = d(i) - c(i) * x(i + 1);
end
end
