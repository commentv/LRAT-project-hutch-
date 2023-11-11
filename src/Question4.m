clear;
d=1000;

f = @(m) 4./(m-2);
g = @(m,d) sqrt(8./((m-4)*d)-1/d^2);

vect_m=8:12:992-2;

dim=1000;
I=eye(dim);
A=I;
repet=100;
mat_rel_err=zeros(1,length(vect_m));

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
    mat_rel_err(1,j)=sqrt(expectation/repet)/answer;
end


m=5:1000;
func1=f(m);
func2=g(m,d);

figure;
ax_1=subplot(1,1,1,'XScale', 'log', 'YScale', 'log');
hold(ax_1,'on');
title(ax_1,'Plot of f_{1}, f_{2} and ARMSE');
xlabel(ax_1,'m');
ylabel(ax_1,'f_{1}, f_{2}, ARMSE');
loglog(ax_1,m,func1,'-b');
loglog(ax_1,m,func2,'-r');
loglog(ax_1,vect_m,mat_rel_err(1,:),'-k');
legend(ax_1,'f_{1}','f_{2}','ARMSE');
legend(ax_1,'Location','northeast')


% Function that computes one iteration of (t-tr(A))^2 %
function E=expectation_one_iterate(dim,m,A,answer,I)
    S=normrnd(0,1,[dim,floor((m+2)/4)]);
    G=randsrc(dim,floor((m-2)/2));
    [Q,~]=qr(A*S,'econ');
    W=(I-Q*Q')*G;
    E=(trace(Q'*A*Q)+2/(m-2)*trace(W'*A*W)-answer)^2;
end
