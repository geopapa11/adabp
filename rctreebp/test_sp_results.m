% PAIRWISE GRAPHS

rand('state',0); randn('state',0); 
Nt=10;  d=3;
Ns=fix(logspace(2,3,Nt)); times=zeros(1,Nt);
Ns=fix(logspace(1,2,Nt)); times=zeros(1,Nt);

for iter=1:100,
for t=1:Nt,

 N=Ns(t);
 A=sparse(N,5*N);
 for i=1:4*N, 
   tmp=1+fix(rand(1,2)*N); A(tmp,i)=1; 
   dims=[d*ones(1,length(unique(tmp))),1]; 
   F{i}=squeeze(rand(dims)); 
 end;

 for i=1:N, A(i,4*N+i)=1; F{4*N+i}=rand(d,1); end;
 At=spanningfg(A);  fkeep=find(sum(At,1)>0); At=At(:,fkeep); F=F(fkeep);
 Ns2(t)=sum(size(At));

 T=rctreeSP(At,F);
 for tests=1:100, 
  Bel=sumproduct(At,F,T.roots);

  for i=1:size(At,1), 
    B2=Bel{i}(:); if (length(B2)==1), B2=ones(d,1)/d; end;
    if (any( (B2(:)-marginal(T,i))./B2(:) > 1e-6 )) fprintf('%d ',i); return; end; 
  end;

  for n=1:20, ind=1+fix(rand*length(fkeep));
   F{ind} = rand(size(F{ind}));
   T=updateFactorSP(T,ind,F{ind});
  end;
  for n=1:20, 
   add=rand > .5; indF=1+fix(rand*length(fkeep)); indV=1+fix(rand*N);
   if (add && length(T.roots)>1)
     F{indF} = rand(d,d); indV=T.roots(1:2);
     if (indV(1)>N), tmp=find(At(:,indV(1)-N)); indV(1)=tmp(1); end;
     if (indV(2)>N), tmp=find(At(:,indV(2)-N)); indV(2)=tmp(1); end;
     At(:,indF)=zeros(N,1); At(indV,indF)=1; 
     if (indV(1)>N), indV(1)=T.Adj{indV(1)}(1); end;
     if (indV(2)>N), indV(2)=T.Adj{indV(2)}(1); end;
     indV=sort(indV);
     [T,Nchanged]=swapFactorSP(T,indF,indV,F{indF});
   else
     F{indF} = rand(d,1);
     At(:,indF)=zeros(N,1); At(indV,indF)=1;
     [T,Nchanged]=swapFactorSP(T,indF,indV,F{indF});
   end;	 
   if (size(At,1)>100), indF,indV, return; end;
  end;
 fprintf('.');
 end;

end;
 fprintf('\n'); 
end;


