% PAIRWISE GRAPHS

rand('state',0); randn('state',0); 
Nt=8;  d=3;
Ns=fix(logspace(2,3,Nt));
Niter=10;
Niter=1;
Nmax=1000;

time1=zeros(Nt,Niter); time2=zeros(Nt,Niter); time3=zeros(Nt,Niter,Nmax); time4=zeros(Nt,Niter,Nmax);
  Nchanged3=zeros(Nt,Niter,Nmax); Nchanged4=zeros(Nt,Niter,Nmax); 

for iter=1:Niter,
for t=1:Nt,
 N=Ns(t);  fprintf('.');
 A=sparse(N,5*N);
 for i=1:4*N, 
   tmp=1+fix(rand(1,2)*N); A(tmp,i)=1; 
   dims=[d*ones(1,length(unique(tmp))),1]; 
   F{i}=squeeze(rand(dims)); 
 end;
 
 for i=1:N, A(i,4*N+i)=1; F{4*N+i}=rand(3,1); end;
% figure(1); imagesc(A);
 At    = spanningfg(A);
 fkeep = find(sum(At,1)>0);
 At    = At(:,fkeep);
 F     = F(fkeep);
 Ns2(t)= sum(size(At));

 % CONSTRUCTION TIME
 tic;
 T=rctreeMP(At,F);
 time1(t,iter)=toc;

 % QUERY TIME
 tic;
 for n=1:Nmax, ind=1+fix(rand*length(fkeep));
   B=maxmarginal(T,ind);
   time2(t,iter,n)=toc; tic;
 end;

 % UPDATE TIME
 tic;
 for n=1:Nmax, ind=1+fix(rand*length(fkeep));
   F{ind} = rand(size(F{ind}));
   T=updateFactorMP(T,ind,F{ind});
   Nchanged3(t,iter,n)=T.changedCount;
   T=clearChanges(T);
   time3(t,iter,n)=toc; tic;
 end;

 % RESTRUCTURE TIME
 tic;
 for n=1:Nmax, 
   add=rand > .5; indF=1+fix(rand*length(fkeep)); indV=1+fix(rand*N);
   if (add && length(T.roots)>1)
     F{indF} = rand(d,d); indV=T.roots(1:2);
     if (indV(1)>N), indV(1)=T.Adj{indV(1)}(1); end;
     if (indV(2)>N), indV(2)=T.Adj{indV(2)}(1); end;
     T=swapFactorMP(T,indF,indV,F{indF});
   else
     F{indF} = rand(d,1);
     T=swapFactorMP(T,indF,indV,F{indF});
   end;	 
   Nchanged4(t,iter,n)=T.changedCount;
   T=clearChanges(T);
   time4(t,iter,n)=toc; tic;
 end;
 
end;
end;

Nchanged3=double(Nchanged3); Nchanged4=double(Nchanged4); 
if (0)

for i=1:Nt,
  X=[];Y=[];S=[];
  for j=0:max(Nchanged3(i,:)),
    ind=find(Nchanged3(i,:)==j);
    Y(j+1)=mean(time3(i,ind),2); S(j+1)=std(time3(i,ind),[],2);
  end; X3=0:max(Nchanged3(i,:)); Y3=Y; L3=Y-S; U3=Y+S;
  X=[];Y=[];S=[];
  for j=0:max(Nchanged4(i,:)),
    ind=find(Nchanged4(i,:)==j);
    Y(j+1)=mean(time4(i,ind),2); S(j+1)=std(time4(i,ind),[],2);
  end; X4=0:max(Nchanged4(i,:)); Y4=Y; L4=Y-S; U4=Y+S;
  figure(1); clf; errorbar(X3,Y3,L3,U3,'b-'); hold on; errorbar(X4,Y4,L4,U4,'r-');
  Y=(X4+1).*(1+log(Ns(i)./X4)); Y=Y*(Y4(end)./Y(end));
  plot(X4, Y, 'g-'); pause(1);
end; 
end;

figure(1);
for i=1:Nt,
  X=0:max(Nchanged4(i,:));  
  Ypeak = mean(time4(i, find(Nchanged4(i,:) > .9*X(end))),2);
  Y=(X+1).*(1+log(Ns(i)./X)); Y=Y*(Ypeak./Y(end));
  subplot(4,2,i);
  plot(Nchanged4(i,:),time4(i,:),'r.',Nchanged3(i,:),time3(i,:),'b.',X,Y,'g-','linewidth',2);
  title(['Ns = ', num2str(Ns(i))]);
  axis tight;
%   pause;
end;    

return;

C1=3*time1(end)/Ns2(end); C2=3*time2(end)/log2(Ns2(end));
figure(2); loglog(Ns2,C1*Ns2,'k-',Ns2,C2*log2(Ns2),'k:',Ns2,time1,'ro-',Ns2,time2,'go-',Ns2,time3,'bo-',Ns2,time4,'co-');
legend('O(n) reference','O(log n) reference','Build','Query','Update','Restructure');
%figure(1); plot(log2(Ns),log2(times),'bo-',log2(Ns2),log2(times),'ro-');

%[time1(1:Nt)./Ns2(1:Nt); time1(1:Nt)./(Ns2(1:Nt).*log2(Ns2(1:Nt)))]
%[time2(1:Nt)./Ns2(1:Nt); time2(1:Nt)./(Ns2(1:Nt).*log2(Ns2(1:Nt)))]
%[time3(1:Nt)./Ns2(1:Nt); time3(1:Nt)./(Ns2(1:Nt).*log2(Ns2(1:Nt)))]


