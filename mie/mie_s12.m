function [S1,S2] = mie_s12(m, x, u)
% Computation of Mie Scattering functions S1 and S2
% for complex refractive index m=m'+im",
% size parameter x=k0*a, and u=cos(scattering angle),
% where k0=vacuum wave number, a=sphere radius;
% s. p. 111-114, Bohren and Huffman (1983) BEWI:TDD122
% C. Mätzler, May 2002
%
% Modified by Shuyu Wei to accept matrix input and output.
% Original code is available at http://omlc.org/software/mie/
% m, x, u must be of the same size
%
% This code is licensed under GNU GPL V2.

nmax=round(2+x+4*x^(1/3));
[an,bn]=mie_abcd(m,x);
an = repmat(an,[1,size(u)]);
bn = repmat(bn,[1,size(u)]);
[pin,tin]=mie_pt(u,nmax);
n=(1:nmax).';
n=(2.*n+1)./(n.*(n+1));
n=repmat(n,[1,size(u)]);
pin=n.*pin;
tin=n.*tin;
S1=squeeze(sum(an.*pin+bn.*tin,1));
S2=squeeze(sum(an.*tin+bn.*pin,1));
