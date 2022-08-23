% Rikesh Sharma 180606 

function z=implicitEuler()
    
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
        phi(n)=phi(n-1)/(1-lambda*deltat);
    end
    subplot(3,2,2);
    t=linspace(0,tmax,tmax/deltat);
    z1=phi;
    plot(t,z1,'Color',[0.8500, 0.3250, 0.0980]);
    title('Implicit Euler for \Deltat=0.1 sec');
    xlabel('t (sec)');
    ylabel('\Phi(t) for \Deltat=0.1 sec');
    
    
    deltat=0.6;
    phi=zeros(2,1);
    phi(1)=1;
    nmax= tmax/deltat;
    for n=2:nmax
        phi(n)=phi(n-1)/(1-lambda*deltat);
    end
    subplot(3,2,3);
    z2=phi;
    t=linspace(0,tmax,tmax/deltat);
    plot(t,z2,'Color',[0.9290, 0.6940, 0.1250]);
    title('Implicit Euler for \Deltat=0.6 sec');
    xlabel('t (sec)');
    ylabel('\Phi(t) for \Deltat=0.6 sec');
    
    
     deltat=2.1;
     phi=zeros(2,1);
     phi(1)=1;
     nmax=tmax/deltat;
     for n=2:nmax
         phi(n)=phi(n-1)/(1-lambda*deltat);
     end
     subplot(3,2,4);
     z3=phi;
     t=linspace(0,tmax,tmax/deltat);
     plot(t,z3,'color',[0.4940, 0.1840, 0.5560]);
     title('Implicit Euler for \Deltat=2.1 sec');
     xlabel('t (sec)');
     ylabel('\Phi(t) for \Deltat=2.1 sec');
    
     
    subplot(3,2,5:6);
    deltat1=0.1; deltat2=0.6; deltat3=2.1;
    t0=linspace(0,tmax,tmax/deltat1);
    t1=linspace(0,tmax,tmax/deltat1);
    t2=linspace(0,tmax,tmax/deltat2);
    t3=linspace(0,tmax,tmax/deltat3);
    plot(t0,z0,t1,z1,t2,z2,t3,z3);
    title('Plots of solution for Exact and Implicit Euler for different value of \Deltat');
    xlabel('t (sec)');
    ylabel('\Phi(t)');
    legend({'\phi(t)=e^{-t}','\Deltat=0.1 sec','\Deltat=0.6 sec','\Deltat=2.1 sec'},'Location','northeast');
     
     %z=[z1(1:20),z2(1:20),z3(1:20)];
end
