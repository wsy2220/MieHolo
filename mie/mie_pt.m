function [p,t] = mie_pt(u,nmax)
% pi_n and tau_n, -1 <= u= cosθ <= 1, n1 integer from 1 to nmax
% angular functions used in Mie Theory
% Bohren and Huffman (1983), p. 94 - 95
% C. Mätzler, May 2002
%
% Modified by Shuyu Wei to accept matrix input and output.
% Original code is available at http://omlc.org/software/mie/
% u is a matrix or vector, nmax is a integer
%
% This code is licensed under GNU GPL V2.

p=zeros([nmax,size(u)]);
t=zeros([nmax,size(u)]);
p(1,:,:)=1;
t(1,:,:)=u;
p(2,:,:)=3.*u;
t(2,:,:)=3.*cos(2.*acos(u));
u = t(1,:,:);
for n1=3:nmax
p1=(2*n1-1)./(n1-1).*p(n1-1,:,:).*u;
p2=n1./(n1-1).*p(n1-2,:,:);
p(n1,:,:)=p1-p2;
t1=n1.*u.*p(n1,:,:);
t2=(n1+1).*p(n1-1,:,:);
t(n1,:,:)=t1-t2;
end;
