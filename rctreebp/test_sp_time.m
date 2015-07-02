% PAIRWISE GRAPHS

rand('state',0); randn('state',0); 
Nt=10;  d=3;
Ns=fix(logspace(2,3.5,Nt)); times=zeros(1,Nt);
time0=times; time1=times; time2=times; time3=times; time4=times;
%Ns=fix(logspace(2,3.6,Nt)); times=zeros(1,Nt);

%rand('state',0); randn('state',0); 
%Nt=1; Ns=100;

%global T
Niter=10;
for iter=1:Niter,
for t=1:Nt,
   if iter>2, t, rand('state'), randn('state'), end;
 N=Ns(t);  fprintf('a');
 Dims = 5.^2+zeros(1,N);
 A=sparse(N,5*N);
 for i=1:4*N, 
   tmp=1+fix(rand(1,2)*N); A(tmp,i)=1; 
   dims=[d*ones(1,length(unique(tmp))),1]; 
   F{i}=squeeze(rand(dims)); 
 end;

 for i=1:N, A(i,4*N+i)=1; F{4*N+i}=rand(d,1); end;
% figure(1); imagesc(A);
 At=spanningfg(A);  fkeep=find(sum(At,1)>0); At=At(:,fkeep); F=F(fkeep);
 Ns2(t)=sum(size(At));

 % CONSTRUCTION TIME
 tic;
 T=rctreeSP(At,F);
 time1(t)=time1(t)+toc;

 fprintf('b'); tic;
 sumproductTree(At,F,T.roots);
 time0(t)=time0(t)+toc/2;

 % QUERY TIME
 fprintf('c'); tic;
 for n=1:100, ind=1+fix(rand*length(fkeep));
   B=marginal(T,ind);
 end;
 time2(t)=time2(t)+toc;

 % UPDATE TIME
 fprintf('d'); tic;
 for n=1:100, ind=1+fix(rand*length(fkeep));
   F{ind} = rand(size(F{ind}));
   T=updateFactorSP(T,ind,F{ind});
 end;
 time3(t)=time3(t)+toc;

 % RESTRUCTURE TIME
 fprintf('e'); tic;
%profile on;
for n=1:100, 
  if (iter>2 && n>0) na=sprintf('break%d.mat',n); rstate=rand('state'); rnstate=randn('state'); save(na); end;
   add=rand > .5; indF=1+fix(rand*length(fkeep)); indV=1+fix(rand*N);
   if (add && length(T.roots)>1)
     F{indF} = rand(d,d); indV=T.roots(1:2);
     if (indV(1)>N), indV(1)=T.Adj{indV(1)}(1); end;
     if (indV(2)>N), indV(2)=T.Adj{indV(2)}(1); end;
     indV=sort(indV);
     [T,Nchanged]=swapFactorSP(T,indF,indV,F{indF});
   else
     F{indF} = rand(d,1);
     [T,Nchanged]=swapFactorSP(T,indF,indV,F{indF});
   end;	 
 end;
%profile report;
 time4(t)=time4(t)+toc/Nchanged;
 
end;
 fprintf('\n'); 
end;
time0=time0/Niter; time1=time1/Niter; time2=time2/Niter; time3=time3/Niter; time4=time4/Niter;
time2=time2/100; time3=time3/100; time4=time4/100;
C1=3*time1(end)/Ns2(end); C2=3*time2(end)/log2(Ns2(end));
%figure(1); clf; loglog(Ns2,C1*Ns2,'k-',Ns2,C2*log2(Ns2),'k:',Ns2,time1,'ro-',Ns2,time2,'go-',Ns2,time3,'bo-',Ns2,time4,'co-');
%legend('O(n) reference','O(log n) reference','Build time','Query time','Update time','Restructure');
figure(1); clf; H=loglog(Ns2,time0,'k-',Ns2,time1,'ro-',Ns2,time2,'go-',Ns2,time3,'bo-',Ns2,time4,'co-');
legend('Sum-product','Build time','Query time','Update time','Restructure');
set(H,'linewidth',2);
%figure(1); plot(log2(Ns),log2(times),'bo-',log2(Ns2),log2(times),'ro-');

%[time1(1:Nt)./Ns2(1:Nt); time1(1:Nt)./(Ns2(1:Nt).*log2(Ns2(1:Nt)))]
%[time2(1:Nt)./Ns2(1:Nt); time2(1:Nt)./(Ns2(1:Nt).*log2(Ns2(1:Nt)))]
%[time3(1:Nt)./Ns2(1:Nt); time3(1:Nt)./(Ns2(1:Nt).*log2(Ns2(1:Nt)))]


