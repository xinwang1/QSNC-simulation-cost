%  Computes the difference between
%  the one-shot QSNC simulation cost and
%  two-shot  average QSNC simulation cost 
%  of a special class of non-commutative bipartite graphs

% Need to install: 
% CVX >= 2.1
% QETLAB >= 0.7

% author: Xin Wang


%   Reference:
%   [1] Xin Wang, Runyao Duan, On the quantum no-signalling assisted 
%       zero-error simulation cost of non-commutative bipartite graphs, 2016 
%   [2] R. Duan, A. Winter, Non-Signalling Assisted Zero-Error Capacity of 
%       Quantum Channels and an Information Theoretic Interpretation of the
%       Lov\'{a}sz Number, arXiv:1409.3426
%   [3] CVX - (http://cvxr.com/cvx/)
%   [4] QETLAB v 0.7 (http://qetlab.com)

function [D] = Simulationcost_difference(r)
%% set dimension 
a=2;b=3;d=a*b;
dim=[a,b];
%% set channel (based on the example in Theorem 9 [1])
a0=sqrt(r);
b0=sqrt(1-a0^2);
v1=[sqrt(1/3) sqrt(1/3) 0 0 0 sqrt(1/3)]';
v2=[0 0 a0 0 b0 0]';
v3=[0 0 0 1 0 0]';
P=v1*v1'+v2*v2'+v3*v3';
%% Compute  one-shot QSNC simulation cost
%Use the SDP of QSNC simulation cost in [2]
cvx_begin sdp quiet
    variable V(d,d) hermitian 
    variable T(b,b) hermitian
    s_cost=trace(T)
    minimize(s_cost)
    subject to
        V >= 0;kron(eye(a),T) - V >= 0;
        PartialTrace(V,2,dim) - eye(a) == 0
        V*(eye(d)-P) == 0;
cvx_end;
sc1=s_cost
%% Compute two-shot average QSNC simulation cost
%set the two-shot dimension
d=d^2;
a=a^2;
b=b^2;
P2=kron(P,P);
dim=[2 3 2 3];
sys=[1,3,2,4];
P2=PermuteSystems(P2,sys,dim);
dim=[a b];

%Use the SDP of QSNC simulation cost in [2]
cvx_begin sdp quiet
    variable V(d,d) hermitian 
    variable T(b,b) hermitian
    s_cost2=trace(T)
    minimize(s_cost2)
    subject to
        V >= 0;kron(eye(a),T) - V >= 0;
        PartialTrace(V,2,dim) - eye(a) == 0
        V*(eye(d)-P2) == 0;
cvx_end
% compute the two-shot average QSNC simulation cost
s_cost2=s_cost2^0.5;

%% Compute the difference
D=s_cost-s_cost2;
end

