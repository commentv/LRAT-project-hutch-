clear;

%Define the parameters%
vect_c=[-0.1,-1,-5];

vect_m=9:9:495;
vect_m(end+1)=498;

dim=1000;
I=eye(dim);
U=create_U(dim);

repet=100;
mat_rel_err=zeros(3,length(vect_m));

for i=1:3
    c=vect_c(i);
    A=U*diag((1:dim).^c)*U';
    answer=trace(A);
    
    % Iteration over m%
    for j=1:length(vect_m)
        j
        m=vect_m(j);
        expectation=0;
        % Iteration over the number of repetition to get an average of
        % (t-tr(A))^2, which would be close to the expectation
        for k=1:repet
            expectation=expectation+expectation_one_iterate(dim,m,A,answer,I);
        end
        % relative error given an m and c %
        mat_rel_err(i,j)=sqrt(expectation/repet)/answer;
    end
end

%Commands to plot%
figure;
ax_1=subplot(1,1,1,'XScale', 'log', 'YScale', 'log');
hold(ax_1,'on');
title(ax_1,'Average relative mean-square error of Hutch++');
xlabel(ax_1,'Number of queries m');
ylabel(ax_1,'Average relative mean-square error \epsilon');
loglog(ax_1,vect_m,mat_rel_err(1,:),'-b');
loglog(ax_1,vect_m,mat_rel_err(2,:),'-r');
loglog(ax_1,vect_m,mat_rel_err(3,:),'-k');
legend(ax_1,'c=0.1','c=1','c=5');
legend(ax_1,'Location','northeast')

% Function that creates a random orthogonal matrix, for the creation of A %
function U=create_U(n)
    rand=normrnd(0,1,[n,n]);
    [U,~]=qr(rand);
end

% Function that computes one iteration of (t-tr(A))^2 %
function E=expectation_one_iterate(dim,m,A,answer,I)
    S=normrnd(0,1,[dim,floor(m/3)]);
    G=normrnd(0,1,[dim,floor(m/3)]);
    [Q,~]=qr(A*S,'econ');
    W=(I-Q*Q')*G;
    E=(trace(Q'*A*Q)+3/m*trace(W'*A*W)-answer)^2;
end