function implicitB()
Re=50;
nu=1/Re;

nx=201;
dx=2/(nx-1);

t=0.47;
dt=0.01;
n=t/dt;

s=(nu*dt)/(dx*dx);

x=-1:dx:1;
u=zeros(201,47);
u(1:101,1)=1;
u(102:201,1)=0;
for i=1:n
    b(1)=1.0; b(47)=1;
    a(47)=0; c(1)=0;
    d(1)=1; d(47)=0;
    
    a(2:nx-1)= -(0.25*(dt/dx)*u(1:nx-2)) - (0.5*s);
    b(2:nx-1)=1+s;
    c(2:nx-1)=(0.25*(dt/dx)*u(3:nx)) - (0.5*s);
    d(2:nx-1)=(0.5*s*u(1:nx-2)) + ((1-s)*u(2:nx-1)) +  (0.5*s*u(3:nx));
    a(nx)=0;
    b(nx)=0;
    c(nx)=0;
    d(nx)=0;
    u(:,i+1)=TDMAx(a,b,c,d,u,nx);
end
for i=1:45
h=plot(x,u(:,i));
hold on
pause(0.1);
delete(h);
end

end

function ansV=TDMAx(a,b,c,d,u,nx)

c(1)=c(1)/b(1);
d(1)=d(1)/b(1);

for i=2:nx
    c(i)=c(i)/(b(i)-a(i)*c(i-1));
    d(i)=(d(i)-a(i)*d(i-1))/(b(i)-a(i)*c(i-1));
end
ansV(nx)=d(nx);
for i=nx-1:-1:1
    ansV(i)=d(i)-u(i+1)*c(i);
end

end
