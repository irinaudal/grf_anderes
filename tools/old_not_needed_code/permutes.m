 function [p,s]=permutes(n)
 %PERMUTES All N! permutations of 1:N + signatures [P,S]=(N)
 % The output P is a matrix of size (N!,N) where each row
 % contains a permutation of the numbers 1:N. The rows are in
 % lexically sorted order.
 %
 % To permute the elements of an arbitrary vector V use
 % V(PERMUTES(LENGTH(V))).

 % PERMUTES(N) is the same as SORTROWS(PERMS(1:N)) but much faster.
 
 % Thanks to Peter J Acklam for several improvements.
 
 %      Copyright (c) 1998 Mike Brookes,  mike.brookes@ic.ac.uk
 %      Version: $Id: permutes.m 713 2011-10-16 14:45:43Z dmb $
 %
 %   VOICEBOX is a MATLAB toolbox for speech processing.
 %   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
 %

 p=1;
 m=1;
 if n>1
   for a=2:n
     q=zeros(a*m,a);
     r=2:a+1;
     ix=1:m;
     for b=1:a
       q(ix,1)=b;
       q(ix,2:a)=r(p);
       r(b)=b;
       ix=ix+m;
     end
     m=m*a;
     p=q;
   end
 end
 if nargout>1 s=1-2*rem(fix((1:m)'/2),2); end