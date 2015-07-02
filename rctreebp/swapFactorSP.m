function [T,Nchanged]=swapFactorSP(T,fid,nbrs,F)

if (nargin < 3), error('Insufficient arguments: swapFactor(T,fid,nbrs,F)'); end;
N=length(T.Dim); NN=length(T.Latest); M=NN-N;
if (fid > M), error('Tree has only %d factors',M); end;
if (any(nbrs > N)), error('Tree has only %d variables',N); end;
fid = uint32(fid+N);  % convert to vertex numbering
nbrs = sort(uint32(nbrs));

% Remove edge (*,fid) from inital graph rep & mark neighbors of vf 
Marked=T.Adj{fid}; 
for i=Marked, modCell(T.Adj,i,sDiff32(T.Adj{i},fid)); end;
% Add edges (nbrs,fid) to inital graph rep & mark nbrs
Marked=sUnion32(Marked,nbrs); 
modCell(T.Adj,fid,nbrs); 
for i=nbrs, modCell(T.Adj,i,sUnion32(T.Adj{i},fid)); end;
Nchanged=length(Marked);
Marked=sUnion32(Marked,fid);
modCell(T.Var,fid,nbrs); modCell(T.F,fid,F);
Reorient=[];

rnd=1;

while (~isempty(Marked)), 
  %fprintf('=== Round %d ===\n',rnd);
  %fprintf('Initial marked: '); for v=Marked, fprintf('%d ',v); end; fprintf('\n');
  for v=Marked,   		% For every marked node, check leaf status & add any other
    flagnew=length(T.Adj{v})<2;  flagold = bitget(T.Roll(v),2);  % affected nodes to list
    %if (v==511), flagnew, flagold, T.Adj{v}, pause; end;
    if (flagnew && ~flagold), Marked=sUnion32(Marked,T.Adj{v}); end; % !!!  
    T.Roll(v)=bitset(T.Roll(v),2,flagnew);
  end;
  Fixup=Marked; MarkNext=Marked;
  %fprintf('Update marked: '); for v=Marked, fprintf('%d ',v); end; fprintf('\n');
  nbrs=Marked;
  for v=Marked,         % For every marked node, re-check the removal decision
    old = ( T.Newer(v)==0 );    % did we remove this node before? & should we now?
    new =~( length(T.Adj{v})>2 || ...  
          (length(T.Adj{v})==2 && (T.Roll(v)==0 || any(T.Roll(T.Adj{v}))) ) || ...
          (length(T.Adj{v})==1 && length(T.Adj{T.Adj{v}})==1 && T.Adj{v}<v )  );
    %fprintf('Node %d: %d -> %d (',v,old,new); for tmp=T.Adj{v}, fprintf('%d ',tmp); end; fprintf(')\n'); pause;
    switch (2*old+new),
    case 0, % (keep,keep) => recompute own link structure in next round
	    Fixup = sUnion32(Fixup,v);
    case {1,3}, % (*,remove) =>  recompute removal, delete newer (if exists), mark nbrs
	    Fixup=sUnion32(Fixup,T.Adj{v}); MarkNext=sUnion32(MarkNext,T.Adj{v}); vv=T.Id(v); 
            T.roots = sDiff32(T.roots,vv); T.Latest(vv)=v;
	      
            remv=T.Newer(v); T.Newer(v)=0; r=rnd;
            while(remv~=0),   %REMOVE v FROM LIST STRUCTURE & FREE MEM
              if (T.Prev(remv)~=0), T.Next(T.Prev(remv))=T.Next(remv); end;
              if (T.Next(remv)~=0), T.Prev(T.Next(remv))=T.Prev(remv); end;
              if (T.First(r)==remv), T.First(r)=T.Next(remv); end;
              T.Next(remv)=T.MemPtr; T.MemPtr=remv; 
              modCell(T.Scars,remv,[]); T.Id(remv)=0; T.Prev(remv)=0; T.Older(remv)=0;
              %fprintf('Freed node %d\n',T.MemPtr);
              remv=T.Newer(remv); r=r+1;
            end;
            nbrs=sUnion32(nbrs,T.Adj{v});

    case 2, % (remove, keep) =>  make & insert newer, update nbrs link structure, mark nbrs
            if (T.MemPtr==0)                        % OUT OF MEMORY?
              %fprintf('Out of memory...\n'); 
              K=T.K; T.K=K+1; NN=length(T.in); KN=K*NN; KKN=KN+NN;	
              T.Next(KKN)=0; T.Prev(KKN)=0;		% If so, extend everything by NN 
              T.Newer(KKN)=0; T.Older(KKN)=0; 
              T.Id(KKN)=0; T.Roll(KKN)=0; 
              T.Adj{KKN}=[];
              T.Next(KN:KKN-1)=KN+1:KKN; T.MemPtr=KN; 
              end;            
            new=T.MemPtr; T.MemPtr=T.Next(T.MemPtr);
            %fprintf('Allocating memory: node %d\n',new);
            T.Next(new)=T.First(rnd+1); T.First(rnd+1)=new;
            if (T.Next(new)~=0), T.Prev(T.Next(new))=new; end;
            T.Latest(T.Id(v))=new; T.Newer(v)=new; T.Older(new)=v; T.Newer(new)=0;
            T.Roll(new)=rand>.5; T.Id(new)=T.Id(v);
	    Fixup=sUnion32(Fixup,sort([v,T.Adj{v}])); MarkNext=sUnion32(MarkNext,T.Adj{v});
            nbrs=sUnion32(nbrs,T.Adj{v});
    end; 
    %fprintf('To Fix: '); for v=Fixup, fprintf('%d (%d) ;',v,T.Newer(v)); end; fprintf('\n');
  end;

  %nbrs=[]; for v=Marked, if (T.Newer(v)==0), nbrs=sUnion32(nbrs,T.Adj{v}); end; end;
  %nbrs=sUnion32(nbrs,Marked);
  for i=nbrs, if (T.Newer(i)~=0), ii=T.Newer(i); modCell(T.Scars,ii,T.Scars{i}); 
    for j=T.Adj{i}, if (T.Newer(j)==0), modCell(T.Scars,ii,sUnion32(sDiff32(T.Scars{ii},T.Scars{j}),j)); end; end;
  end; end;

  for v=Marked, if (T.Newer(v)==0),
    %% deleting the node -- get our scars, remove from circulation, add self
    vv=T.Id(v); cfs=T.Scars{v}; nbrs=T.Adj{v};
    vars=T.Var{vv}; for i=T.Id(cfs), 
      if (i<1 || i~=real(round(i))), v,cfs,T.Id(cfs),nbrs, pause; end;
      vars=sUnion32(vars,T.CFVar{i}); 
    end; dims=T.Dim(vars);
    CF = multCF(1,vars,dims,T.F{vv},T.Var{vv});            
    for i=T.Id(cfs), CF = multCF(CF,vars,dims,T.CF{i},T.CFVar{i}); end;
    switch (length(nbrs))
    case 0, CF=sum(CF(:)); s1=[]; T.roots=sUnion32(T.roots,vv); keep=[]; T.in(vv)=0; T.out(vv)=0;
    case 1, T.in(vv)=T.Id(T.Adj{v}); T.out(vv)=0; keep=T.Var{T.Id(nbrs)};
    case 2, T.in(vv)=T.Id(nbrs(1)); T.out(vv)=T.Id(nbrs(2)); 
            Reorient=sUnion32(Reorient,vv);
            keep=sUnion32(T.Var{T.Id(nbrs(1))},T.Var{T.Id(nbrs(2))});
    end;
    hits=sMember32(vars,keep); i1=find(hits); s1=vars(i1); i2=find(~hits);
    if (~isscalar(vars)), CF=permute(CF,[i1,i2]); end;
    for i=length(vars):-1:length(s1)+1, CF=sum(CF,i); end;
    modCell(T.CF,vv,CF); modCell(T.CFVar,vv,s1);
  end; end;

  Fixup=sDiff32(sort(T.Newer(Fixup)),uint32(0)); MarkNext=sDiff32(sort(T.Newer(MarkNext)),uint32(0)); rnd=rnd+1;
  for ptr=Fixup,   % Fix up the adjacency structure of the new round
    old=T.Older(ptr); AdjO=T.Adj{old}; AdjN=T.Newer(AdjO);
    for i=find(AdjN==0),
      tmp = T.Newer(T.Adj{AdjO(i)}); 
      if (length(tmp)>0 && tmp(1)~=ptr), AdjN(i)=tmp(1); elseif (length(tmp)>1), AdjN(i)=tmp(2); end;
    end;
    AdjN=AdjN(AdjN~=0); modCell(T.Adj,ptr,AdjN);
  end;
  Marked = MarkNext;
  
end; 

% Orientation now done during query operation...
%%%%%%
%todo=zeros(1,NN); todo(Reorient)=1; 
%while (~isempty(Reorient)), v=Reorient(end);
% if todo(v),
%  u=T.in(v); w=T.out(v);
%  if (todo(u)), Reorient=[Reorient,u]; continue; end;
%  if (todo(w)), Reorient=[Reorient,w]; continue; end;
%  if ((T.in(u)==w)||(T.out(w)==u)), T.in(v)=w; T.out(v)=u; else T.in(v)=u; T.out(v)=w; end;
%  todo(v)=0;
% end; 
% Reorient=Reorient(1:end-1);
%end;
%

