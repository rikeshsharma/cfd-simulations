function exactSol()
Re=10;
nu=1/Re;
xmax=1;
dx=0.01;
t=0;
x=-xmax:dx:xmax;
N=2*xmax/dx+1;
u=zeros(N,100);
for j=1:100
    t=t+0.01*j;
    for i = 1:N
    
        f = @(z)  ((x(i)-z)/t).*exp(-0.5*Re*(z+0.5*((x(i)-z).^2)/t)) ;
        g = @(z)  ((x(i)-z)/t).*exp(-0.5*Re*(0.5*((x(i)-z).^2)/t)) ;
        h = @(z)  exp(-0.5*Re*(z+0.5*((x(i)-z).^2)/t));
        l = @(z)  exp(-0.5*Re*(0.5*((x(i)-z).^2)/t));
    
        fi = integral(f,-inf,0);
        gi = integral(g,0,inf);
        hi = integral(h,-inf,0);
        li = integral(l,0,inf);
    
        u(i,j)=(fi+gi)/(hi+li);
    end
end


xlabel(' X from (-1 to 1) ');
ylabel('U for x at time t=0.47 sec');
title('RIKESH SHARMA 180606','U for x at time t=0.47 sec');

for i=1:30
    h=plot(x,u(1:N,i),'r');
    hold on
    pause(0.2);
    delete(h);
end
    h=plot(x,u(1:N,30),'r');
end