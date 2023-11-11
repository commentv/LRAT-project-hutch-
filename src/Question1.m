clear;

%Define the parameters%
dim=1000;
rep=100;
a=1:dim;
l=9:9:495;
l(end+1)=498;
Q=genort(dim);
L=zeros([3,length(l)]);
c1=[0.1,1,5];

%Iteration over indices of c1, c1(i) indicates value of c
for k=1:3
    c=c1(k);
    A=diag(a.^(-c));
    W=Q*A*Q.';
    %s - real value of trace
    s=sum(a.^(-c));
    S=zeros([1,length(l)]);
    %Iteration over all m=l(i) values
    for i=1:length(l)
        i
        %Iteration over the realization number to calculate
        %the estimation of the mean of square differences of
        %the l(i)-Hutchinson estimator and real trace s.
        for j=1:rep
            S(i)=S(i)+approx(W,l(i),s,dim);
        end
        S(i)=S(i)/rep;
    end
    %Average relative errors are put into a matrix L for plotting
    L(k,:)=sqrt(S)/s;
end

%Plotting the average relative error
figure;
ax_1=subplot(1,1,1,'XScale', 'log', 'YScale', 'log');
hold(ax_1,'on');
title(ax_1,'Average relative mean-square error of Monte-Carlo estimator of the trace');
xlabel(ax_1,'Number of queries m');
ylabel(ax_1,'Average relative mean-square error \epsilon');
loglog(ax_1,l,L(1,:),'-b');
loglog(ax_1,l,L(2,:),'-r');
loglog(ax_1,l,L(3,:),'-m');
legend(ax_1,'c=0.1','c=1','c=5');
legend(ax_1,'Location','northeast')

%Function generating an orthogonal matrix as the
%Jordan normal form of a symmetric matrix
function Q=genort(dim)
    Q1=normrnd(0,1,[dim,dim]);
    [Q,~]=eig(Q1+Q1.');
end

%Function calculating the square difference between
%the m-Hutchinson estimator for the current
%realization and the real trace value 
function ans=approx(W,m,s,dim)
    sum=0;
    for i=1:m
        w=normrnd(0,1,[1,dim]);
        sum=sum+(w*W*w.'-s);
    end
    ans=(sum/m)^2;
end
    
    