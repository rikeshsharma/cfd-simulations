% Rikesh Sharma 180606 

function z=explicitEuler()
    %Exact solution
    w2=1;
    w=1;
    phi0=1;
    deltat=0.1;
    tmax=15; 
    t=linspace(0,tmax,tmax/deltat);
    phi=phi0*cos(w*t);
    z0=phi;
    subplot(3,2,1);
    plot(t,z0);
    title('Exact solution');
    xlabel('t (sec)');
    ylabel('\Phi(t)=\phi_ocos(\omegat)');
    hold on
     
    % for delat t=0.1
    phi=zeros(2,1);
    phid=zeros(2,1);
    phi(1)=1;
    phid(1)=0;
    nmax=tmax/deltat;
    for n=2:nmax
        phid(n)=phid(n-1)-w2*phi(n-1)*deltat;
        phi(n)=phi(n-1)+phid(n-1)*deltat;
    end
    subplot(3,2,2);
    z1=phi;
    t=linspace(0,tmax,tmax/deltat);
    plot(t,z1,'Color',[0.8500, 0.3250, 0.0980]);
    title('Explicit Euler for \Deltat=0.1 sec');
    xlabel('t (sec)');
    ylabel('\Phi(t) for \Deltat=0.1 sec');
   
    % for delat t=0.6
    deltat=0.6;
    phi=zeros(2,1);
    phid=zeros(2,1);
    phi(1)=1;
    phid(1)=0;
    nmax=tmax/deltat;
    for n=2:nmax
        phid(n)=phid(n-1)-w2*phi(n-1)*deltat;
        phi(n)=phi(n-1)+phid(n-1)*deltat;
    end
    subplot(3,2,3);
    z2=phi;
    t=linspace(0,tmax,tmax/deltat);
    plot(t,z2,'Color',[0.9290, 0.6940, 0.1250]);
    title('Explicit Euler for \Deltat=0.6 sec');
    xlabel('t (sec)');
    ylabel('\Phi(t) for \Deltat=0.6 sec');
   
    % for delat t=2.1
    deltat=2.1;
    phi=zeros(2,1);
    phid=zeros(2,1);
    phi(1)=1;
    phid(1)=0;
    nmax=tmax/deltat;
    for n=2:nmax
        phid(n)=phid(n-1)-w2*phi(n-1)*deltat;
        phi(n)=phi(n-1)+phid(n-1)*deltat;
    end
    subplot(3,2,4);
    z3=phi;
    t=linspace(0,tmax,tmax/deltat);
    plot(t,z3,'color',[0.4940, 0.1840, 0.5560]);
    title('Explicit Euler for \Deltat=2.1 sec');
    xlabel('t (sec)');
    ylabel('\Phi(t) for \Deltat=2.1 sec');
    
    %all plots in one graph
    subplot(3,2,5:6);
    deltat1=0.1; deltat2=0.6; deltat3=2.1;
    t0=linspace(0,tmax,tmax/deltat1);
    t1=linspace(0,tmax,tmax/deltat1);
    t2=linspace(0,tmax,tmax/deltat2);
    t3=linspace(0,tmax,tmax/deltat3);
    plot(t0,z0,t1,z1,t2,z2,t3,z3);
    title('Plots of solution for Exact and Explicit Euler for different value of \Deltat');
    xlabel('t (sec)');
    ylabel('\Phi(t)');
    legend({'\Phi(t)=\phi_ocos(\omegat)','\Deltat=0.1 sec','\Deltat=0.6 sec','\Deltat=2.1 sec'},'Location','northwest');
  % z=[z1(1:20),z2(1:20),z3(1:20)];
   
end