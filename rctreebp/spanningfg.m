

function [AA]=spanningfg(A)

[N,M]=size(A);
AA=sparse(N,M);
uV=zeros(1,N); uF=zeros(1,M); 
for i=1:M
  fnd=find(uF<2);
  if (~isempty(fnd))
   ind=fnd(1+fix(rand(1)*length(fnd))); uF(ind)=3;
   AA(:,ind)=A(:,ind);	% pick random factor, add var edges
   fnd=find(A(:,ind))'; uV(fnd)=1; % mark vars as used,
   for j=fnd, nbr=find(A(j,:)); uF(nbr)=uF(nbr)+1; end;
  end;
  %figure(2); imagesc(AA); figure(3); plot(uF); pause;
end;


