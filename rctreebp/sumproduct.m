% [Bel,Msgs] = sumproduct(Adj,Fac,roots)
%   a quick & dirty BP implementation for factor graphs.  For N vars & M 
%   factors, Adj is a NxM matrix of 0/1 and Fac is 1xM cell array of factors
%   (matrices) of appropriate dimension.  Runs BP slightly inefficiently:
%   computes a breadth-first traversal of the tree from the list "roots",
%   then messages from leaves to root & back, but each message is actually computed
%   twice (once on each pass) rather than the minimal once.
%
function [Bel,Msgs] = sumproduct(Adj,Fac,roots) 

[N,M]=size(Adj);
order=zeros(1,N+M); par=zeros(1,N+M);
order(1:length(roots))=roots; j=length(roots)+1;
for i=1:N+M,
  %if (order(i)==0) order=order(1:i-1); break; end;
  if (order(i)>N), nbrs=find(Adj(:,order(i)-N))';
  else             nbrs=find(Adj(order(i),:))+N;
  end;
  nbrs = setdiff_sorted(nbrs,par(i)); Ln=length(nbrs);
  order(j:j+Ln-1)=nbrs; par(j:j+Ln-1)=order(i); j=j+Ln;
  %fprintf('%d ',i);
end;
order = [order(end:-1:1),order];

[N,M] = size(Adj); Dims=zeros(1,N); fnbrs=cell(1,N); vnbrs=cell(1,M);
for i=1:N, fnbrs{i}=find(Adj(i,:)); end; 	% convert to more useful
for i=1:M, vnbrs{i}=find(Adj(:,i))'; 
  tmp=size(Fac{i}); if (tmp(end)==1) tmp=tmp(1:end-1); end;
  if (length(vnbrs{i})>0), Dims(vnbrs{i})=tmp; end;
end; 						%   neighborhood structure

Msgs = cell(N+M); Bel = cell(1,N);
for j=1:M,
  for ii=1:length(vnbrs{j}), i=vnbrs{j}(ii);
    Msgs{i,N+j} = ones(1,size(Fac{j},ii));
    Msgs{N+j,i} = Msgs{i,N+j};
end; end;

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for it=order
	
if (it <= N), i=it;
  Bel{i} = 1;			%  and beliefs at each variable node
  for j=N+fnbrs{i},
    Bel{i} = Bel{i} .* Msgs{j,i}; Bel{i}=Bel{i}./sum(Bel{i});
  end;
  for j=N+fnbrs{i}
    Msgs{i,j} = Bel{i} ./ Msgs{j,i}; Msgs{i,j}=Msgs{i,j}./sum(Msgs{i,j});
  end; 
end; 

if (it > N), j=it-N;
 m=[1]; for ii=1:length(vnbrs{j}), i=vnbrs{j}(ii); 
   m2=Msgs{i,N+j};
   dims=size(m2,2); 	% !!! ugly.  clean up with multCF function...
   n=size(m,1);
   m = repmat(m,[1,dims]) .* repmat(m2,[n,1]);
   m = reshape(m,[n*dims,1]);
 end; 
 if (length(vnbrs{j})>1), m = reshape(m,Dims(vnbrs{j})); end;

 for ii=1:length(vnbrs{j}), i=vnbrs{j}(ii);
   %(2) convert f(x,y,z) -> f(y ; [x z])
   f = Fac{j}; ndim=length(size(f)); f = permute(f,[ii,1:ii-1,ii+1:ndim]);
   f = reshape(f,size(f,1),numel(f)./size(f,1));
   %(3) convert m(x), m(z) -> m(xz)
   [tmp,ind]=max(Msgs{i,N+j});				% compute m(x)*m(y=ind)*m(z) \propto m(x)*m(z)
   mi = permute(m,[ii,1:ii-1,ii+1:ndim]); mi=mi(ind,:);
   tmp = f*mi'; Msgs{N+j,i} = tmp'/sum(tmp);		% then convolve it with f(y ; [x z])
 end;
end; 

end;
