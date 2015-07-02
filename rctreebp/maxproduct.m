% [Bel,Msgs] = maxproduct(Adj,Fac,roots)
%   a quick & dirty BP implementation for factor graphs.  For N vars & M 
%   factors, Adj is a NxM matrix of 0/1 and Fac is 1xM cell array of factors
%   (matrices) of appropriate dimension.  Runs leaf-to-root & back BP
%   , from initial messages Msgs (uniform if not specified).
%   (note: beliefs are computed using "old" messages, i.e. iter-1 iterations)
%
function [Bel,Msgs] = maxproduct(Adj,Fac,roots) 

%%%%%%%%%%%%%%%%% FIND ORDER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[N,M]=size(Adj);				% find # vars, # factors
order=zeros(1,N+M); par=zeros(1,N+M);		% find an order for processing
order(1:length(roots))=roots; j=length(roots)+1;
for i=1:N+M,
  if (order(i)>N), nbrs=find(Adj(:,order(i)-N))';
  else             nbrs=find(Adj(order(i),:))+N;
  end;
  nbrs = setdiff_sorted(nbrs,par(i)); Ln=length(nbrs);
  order(j:j+Ln-1)=nbrs; par(j:j+Ln-1)=order(i); j=j+Ln;
end;
order = [order(end:-1:1),order];		% leaves first, then roots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%% CONVERT STRUCTURE & ALLOCATE MESSAGES %%%%%%%%%%%%%%%%%%%%%%
[N,M] = size(Adj); Dims=zeros(1,N); fnbrs=cell(1,N); vnbrs=cell(1,M);
for i=1:N, fnbrs{i}=find(Adj(i,:)); end; 	% convert to more useful
for i=1:M, vnbrs{i}=find(Adj(:,i))'; 
  tmp=size(Fac{i}); if (tmp(end)==1) tmp=tmp(1:end-1); end;
  Dims(vnbrs{i})=tmp; 
end; 						%   neighborhood structure

Msgs = cell(N+M); Bel = cell(1,N);
for j=1:M,
  for ii=1:length(vnbrs{j}), i=vnbrs{j}(ii);
    Msgs{i,N+j} = ones(1,size(Fac{j},ii));
    Msgs{N+j,i} = Msgs{i,N+j};
end; end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

for it=order
	
if (it <= N), i=it;		%%% VARIABLE NODE %%%
  Bel{i} = 1;			% Compute beliefs from incoming msgs 
  for j=N+fnbrs{i},		%
    Bel{i} = Bel{i} .* Msgs{j,i}; Bel{i}=Bel{i}./sum(Bel{i});
  end;				% Then divide to find outgoing msgs
  for j=N+fnbrs{i}
    Msgs{i,j} = Bel{i} ./ Msgs{j,i}; Msgs{i,j}=Msgs{i,j}./sum(Msgs{i,j});
  end; 
end; 

if (it > N), j=it-N;		%%% FACTOR NODE %%%
 m=[1]; for ii=1:length(vnbrs{j}), i=vnbrs{j}(ii); 
   m2=Msgs{i,N+j};
   dims=size(m2,2); 		% !!! ugly.  clean up with multCF function...
   n=size(m,1);
   m = repmat(m,[1,dims]) .* repmat(m2,[n,1]);
   m = reshape(m,[n*dims,1]);
 end; 
 if (length(vnbrs{j})>1), m = reshape(m,Dims(vnbrs{j})); end;

 for ii=1:length(vnbrs{j}), i=vnbrs{j}(ii);
   f = Fac{j}; ndim=length(size(f)); 		% convert factor f(x,y,z) -> f(y ; [x z])
   f = permute(f,[ii,1:ii-1,ii+1:ndim]);	%   suitable for msg to "y"
   f = reshape(f,size(f,1),numel(f)./size(f,1));% 
   mi = permute(m,[ii,1:ii-1,ii+1:ndim]);       % convert m(x)*m(y)*m(z) -> m(y; xz)
   [tmp,ind]=max(Msgs{i,N+j}); mi=mi(ind,:);	%   then m(xz) \propto m(y=ind; xz)

%   tmp = f*mi'; 				% sum-product: convolve with f
   tmp = max( f.*repmat(mi,[Dims(i),1]) ,[],2); % max-product: maximize over xz

   Msgs{N+j,i} = tmp'/sum(tmp);
 end;
end; 

end;
