% Rikesh Sharma 180606 

function plotResult()
    
    %Exact Solution
    subplot(4,2,1);
    x=0:0.01:1;
    y_exact=exp(x);
    error=y_exact-y_exact;
    plot(x,error,'linewidth',1);
    title('No Error exact solution');
    xlabel('X');
    ylabel('zero error');
    
    hold on
    
    %for N=1
    subplot(4,2,2);
    A=1/2;
    b=1;
    a=A\b;
    x_vec=x;
    y1=1+a*x_vec;
    error=y_exact-y1;
    plot(x,error,'-rs','MarkerIndices',1:10:length(x),'MarkerEdgeColor',[0, 0.4470, 0.7410]	);
    title('Error for N=1');
    xlabel('X');
    ylabel('Yexact-Y');
    
    disp(a);
    
    %for N=2
    subplot(4,2,3);
    A=[1/2, 2/3; 
        1/6, 5/12];
    b=[1;1/2];
    a=A\b;
    x_vec=[x;x.*x];
    y2=1+a'*x_vec;
    error=y_exact-y2;
    plot(x,error,'-gs','MarkerIndices',1:10:length(x),'MarkerEdgeColor',[0, 0.4470, 0.7410]	);
    title('Error for N=2');
    xlabel('X');
    ylabel('Yexact-Y');
    
    disp(a);
    
    %for N=4
    subplot(4,2,4);
    A = [ 1/2,  2/3,  3/4,   4/5; 
          1/6,  5/12, 11/20, 19/30;
          1/12, 3/10, 13/30, 11/21;
          1/20, 7/30, 5/14,  25/56];
    b = [1; 1/2; 1/3; 1/4];
    a = A\b;
    x_vec = [x; x.*x; x.*x.*x; x.*x.*x.*x];
    y3 = 1+a'*x_vec;
    error=y_exact-y3;
    plot(x,error,'-bs','MarkerIndices',1:10:length(x),'MarkerEdgeColor',[0, 0.4470, 0.7410]	);
    title('Error for N=4');
    xlabel('X');
    ylabel('Yexact-Y');
    
    disp(a);
    
    %for N=5
    subplot(4,2,5);
    A = [ 1/2,  2/3,  3/4,   4/5,   5/6; 
          1/6,  5/12, 11/20, 19/30, 29/42;
          1/12, 3/10, 13/30, 11/21, 33/56;
          1/20, 7/30, 5/14,  25/56, 37/72;
          1/30, 4/21, 17/56, 7/18,  41/90];
    b = [1; 1/2; 1/3; 1/4; 1/5];
    a = A\b;
    x_vec = [x; x.*x; x.*x.*x; x.*x.*x.*x; x.*x.*x.*x.*x];
    y4 = 1+a'*x_vec;
    error=y_exact-y4;
    plot(x,error,'-cs','MarkerIndices',1:10:length(x),'MarkerEdgeColor',[0, 0.4470, 0.7410]	);
    title('Error for N=5');
    xlabel('X');
    ylabel('Yexact-Y');
    
    disp(a);
    
    %for N=7
    subplot(4,2,6);
    A = [ 1/2,  2/3,  3/4,   4/5,   5/6,    6/7,    7/8; 
          1/6,  5/12, 11/20, 19/30, 29/42,  41/56,  55/72;
          1/12, 3/10, 13/30, 11/21, 33/56,  23/36,  61/90;
          1/20, 7/30, 5/14,  25/56, 37/72,  17/30,  67/110;
          1/30, 4/21, 17/56, 7/18,  41/90,  28/55,  73/132;
          1/42, 9/56, 19/72, 31/90, 9/22,   61/132, 79/156;
          1/56, 5/36, 7/30,  17/55, 49/132, 11/26,  85/182];
    b = [1; 1/2; 1/3; 1/4; 1/5; 1/6; 1/7];
    a = A\b;
    x_vec = [x; x.*x; x.*x.*x; x.*x.*x.*x; x.*x.*x.*x.*x; x.*x.*x.*x.*x.*x; x.*x.*x.*x.*x.*x.*x];
    y5 = 1+a'*x_vec;
    error=y_exact-y5;
    plot(x,error,'-ms','MarkerIndices',1:10:length(x),'MarkerEdgeColor',[0, 0.4470, 0.7410]	);
    title('Error for N=7');
    xlabel('X');
    ylabel('Yexact-Y');
    
    disp(a);
    
    %plot of all functions
    subplot(4,2,7:8)
    plot(x,y_exact,x,y1,x,y2,x,y3,x,y4,x,y5);
    title('Plot of solution of dy/dx=y for exact and by FEM');
    xlabel('X');
    ylabel('Y');
    legend({'Exact solution','N=1','N=2','N=4','N=5','N=7'},'Location','northwest');
    
    
end