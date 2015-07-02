% T=updateFactor(T,i,F) : change the ith factor in (the factor graph of) T
%   to a new matrix, F.  F should be the same size as the old factor.
%
function [T,t_upd_map]=updateFactorMP(T,ind,F)

u=0; v=ind+length(T.Dim);     	% convert factor to vertex number (|dim|=N)
modCell(T.F,v,F); %T.F{v}=F;			% save the factor itself
updated=[];
while (v~=0),                 	% from v up to root, update the cluster f'ns

  updated=[updated,v];

  %fprintf('Updating node %d\n',v);
  cfs = T.Id(T.Scars{T.Latest(v)});  % get v's child cluster functions
  vars=T.Var{v};    		     % get variables and their dimensions...
  for i=1:length(cfs), vars = sUnion32(vars, T.CFVar{cfs(i)}); end;
  dims=T.Dim(vars);		% save the dimension of involved variables

  % compute the resulting cluster function (by marginalizing)
  CF = multCF(1,vars,dims,T.F{v},T.Var{v});
  for i=1:length(cfs)
    CF = multCF(CF,vars,dims,T.CF{cfs(i)},T.CFVar{cfs(i)});
  end;
  keep=[]; 
  if (T.in(v)), keep=sUnion32(keep,T.Var{T.in(v)}); end;
  if (T.out(v)), keep=sUnion32(keep,T.Var{T.out(v)}); end;
  hits=sMember32(vars,keep); i1=find(hits); s1=vars(i1); i2=find(~hits);
  %[s1,i1]=intersect(vars,keep); [s2,i2]=setdiff(vars,keep);
  if (~isscalar(vars)) CF=permute(CF,[i1,i2]); end;
  for i=length(vars):-1:length(s1)+1, CF=max(CF,[],i); end;
  modCell(T.CF,v,CF); modCell(T.CFVar,v,s1);	% save the new cluster function 
	
  u=v;  w1=T.in(u); w2=T.out(u);	% find the next parent node
  if ((w2==0)||(T.in(w1)==w2)||(T.out(w1)==w2)) v=w1; else v=w2; end;
end;

t_s_map = tic;
for i=updated(end:-1:1), searchDown(T,i); end;	% check for changes along route
t_upd_map = toc(t_s_map);


