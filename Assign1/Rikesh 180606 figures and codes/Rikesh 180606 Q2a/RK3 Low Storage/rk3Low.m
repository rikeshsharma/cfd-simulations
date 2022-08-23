% Rikesh Sharma 180606 

function z=rk3Low()

    deltat=0.1;
    lambda=-1;
    tmax=10; 
    t=linspace(0,tmax,tmax/deltat);
    phi=exp(-t);
    z0=phi;
    subplot(3,2,1);
    plot(t,z0);
    title('Exact solution');
    xlabel('t (sec)');
    ylabel('\Phi(t)=e^{-t}');
    hold on
    
    phi=zeros(2,1); 
    phi(1)=1;
    nmax=tmax/deltat;
    for n=2:nmax
        phi_t=phi(n-1);
        k1=lambda*phi_t;
        phi_t=phi_t+deltat*k1/3;
        k1=-5*k1/9 + lambda*phi_t;
        phi_t=phi_t+15*deltat*k1/16;
        k1=-153*k1/128 +lambda*phi_t;
        phi(n)= phi_t + 8*deltat*k1/15;
    end
    subplot(3,2,2);
    z1=phi;
    t=linspace(0,tmax,tmax/deltat);
    plot(t,z1,'Color',[0.8500, 0.3250, 0.0980]);
    title('RK3 Low for \Deltat=0.1');
    xlabel('t (sec)');
    ylabel('\Phi(t) for \Deltat=0.1');
    
    
    
    phi=zeros(2,1);
    deltat = 0.6;
    phi(1)=1;
    nmax=tmax/deltat;
    
    for n=2:nmax
        phi_t=phi(n-1);
        k1=lambda*phi_t;
        phi_t=phi_t+deltat*k1/3;
        k1=-5*k1/9 + lambda*phi_t;
        phi_t=phi_t+15*deltat*k1/16;
        k1=-153*k1/128 +lambda*phi_t;
        phi(n)= phi_t + 8*deltat*k1/15;
    end
    subplot(3,2,3);
    z2=phi;
    t=linspace(0,tmax,tmax/deltat);
    plot(t,z2,'Color',[0.9290, 0.6940, 0.1250]);
    title('RK3 Low for \Deltat=0.6');
    xlabel('t (sec)');
    ylabel('\Phi(t) for \Deltat=0.6');
    
    
    phi=zeros(2,1);
    deltat = 2.1;
    phi(1)=1;
    nmax=tmax/deltat;
    
    for n=2:nmax
        phi_t=phi(n-1);
        k1=lambda*phi_t;
        phi_t=phi_t+deltat*k1/3;
        k1=-5*k1/9 + lambda*phi_t;
        phi_t=phi_t+15*deltat*k1/16;
        k1=-153*k1/128 +lambda*phi_t;
        phi(n)= phi_t + 8*deltat*k1/15;
    end
    subplot(3,2,4);
    z3=phi;
    t=linspace(0,tmax,tmax/deltat);
    plot(t,z3,'color',[0.4940, 0.1840, 0.5560]);
    title('RK3 Low for \Deltat=2.1');
    xlabel('t (sec)');
    ylabel('\Phi(t) for \Deltat=2.1');
    
    subplot(3,2,5:6);
    deltat1=0.1; deltat2=0.6; deltat3=2.1;
    t0=linspace(0,tmax,tmax/deltat1);
    t1=linspace(0,tmax,tmax/deltat1);
    t2=linspace(0,tmax,tmax/deltat2);
    t3=linspace(0,tmax,tmax/deltat3);
    plot(t0,z0,t1,z1,t2,z2,t3,z3);
    title('Plots of solution for Exact and RK3 Low for different value of \Deltat');
    xlabel('t (sec)');
    ylabel('\Phi(t)');
    legend({'\phi(t)=e^{-t}','\Deltat=0.1','\Deltat=0.6','\Deltat=2.1'},'Location','northeast');
    %z=[z1(1:20),z2(1:20),z3(1:20)];
end