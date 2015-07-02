% updateMAP(T,v) -- compute the MAP config of v conditioned on the MAP config above v.
%
%
%
function vars=updateMAP(T,v)
  %fprintf('Computing max-marginal %d\n',v);
  orient(T,v,0);	% reorient v without recursion (parents oriented)
  % get req'd messages:
  cVars = [];
  if (T.out(v)~=0) cVars=[cVars,T.Var{T.out(v)}]; end;	% find vars already optimized at
  if (T.in(v)~=0) cVars=[cVars,T.Var{T.in(v)}]; end;	%   parent nodes of T
  cVars = sort(cVars);					% sorted order
  cVals = T.map(cVars);					% get map config of those vars
 
  vars = T.Var{v}; cfs=T.Id(T.Scars{T.Latest(v)});	% get info from below v 
  for i=cfs, vars=sUnion32(vars, T.CFVar{i}); end; dims=T.Dim(vars);

  % compute product of all cluster f'ns, incoming msgs, etc:
  CF = multCF(1,vars,dims,T.F{v},T.Var{v});
  for i=1:length(cfs)
    CF = multCF(CF,vars,dims,T.CF{cfs(i)},T.CFVar{cfs(i)});
  end;

  %fprintf('\n Conditioned max-marg:\n'); vars,cVars, pause;

  varsSave=vars;

  hits=sMember32(vars,cVars); i1=find(hits); i2=find(~hits);
  if (~isscalar(vars)) CF=permute(CF,[i2,i1]); end;
  if (isempty(i2))					% if all variables are set by parents
    vals=[]; vars=[];					% nothing to set here
  else 
    szCF=size(CF); tmp=prod(szCF);			%
    for i=length(i1):-1:1				% otherwise, for each variable set by parents 
      CF=reshape(CF,[tmp/szCF(end),szCF(end)]);		%  reshape so it's the second & last dimension
      CF=CF(:, T.map(vars(i1(i))));			%  condition on its value
      szCF=szCF(1:end-1); tmp=prod(szCF);		%  and remove it from the resizing vector
    end;
    if (length(szCF)==1) szCF=[szCF,1]; end;		% put it back into the shape of the
    CF=reshape(CF,szCF);				%  remaining variables

    vars=vars(i2);					% get only the newly optimized vars
    [mx,indx] = max(CF(:));				% find their maximizing value
    vals=ind2sub2(T.Dim(vars),indx);			% and convert into an index vector

    for i=1:length(vars), ii=uint32(vars(i));		% check to see if any of these are different
      if (T.map(ii)~=vals(i)) 				%  from the saved MAP values
        modU32Array(T.map,ii,uint32(vals(i)));		% if so store the new value & mark them as changed
        modU32Array(T.diff,ii,uint32(1));		%   T.diff(ii)=1       -- mark as changed
        modU32Array(T.changedCount,uint32(1),T.changedCount+1);	%   T.changedCount++   -- increase changed list
        modU32Array(T.changed,T.changedCount,ii);	%   T.changed(cnt)=ii  -- add to changed list
      end;
    end;


  end;
  vars = varsSave;
%  vars=T.Var{v}; vals=T.map(vars);			% return the variables & their map configuration
















%  
% Compute (downward) RC-tree message m_{u->v}
%
function [msg,var]=maxmessageC(T,u,v)
  u=uint32(u); v=uint32(v);
  if (u==0) msg=1; var=[]; return; end;
  %fprintf('Computing max-product conditional message %d -> %d\n',u,v);
  % % follow recursion to get incoming message:
  % if (T.in(v)==u) [msg,var]=maxmessage(T,T.in(u),u);
  % else [msg,var]=maxmessage(T,T.out(u),u);          
  % end;
  tmp=[0,v]; 
  while (tmp(2)~=u), 
    tmp(1)=tmp(2);  uu=T.in(tmp(1)); ww=T.out(tmp(1));
    if ((T.in(uu)==ww)||(T.out(uu)==ww)) tmp(2)=uu; else tmp(2)=ww; end;
  end; chEx=tmp(1);
  % get list of all involved variables:
  vars=T.Var{u}; %vars=sUnion32(var,T.Var{u}); 
  cfs=sDiff32( sort(T.Id(T.Scars{T.Latest(u)})) , chEx );
  for i=cfs, vars=sUnion32(vars, T.CFVar{i}); end;
  dims=T.Dim(vars);
  % compute product of all cluster f'ns, etc.:
  CF=multCF(CF,vars,dims,T.F{u},T.Var{u});
  for i=1:length(cfs)
    CF = multCF(CF,vars,dims,T.CF{cfs(i)},T.CFVar{cfs(i)});
  end;
  keep = T.CFVar{v};   		% maximize over "other" variables
  hits=sMember32(vars,keep); i1=find(hits); s1=vars(i1); i2=find(~hits);
  %[s1,i1]=intersect(vars,keep); [s2,i2]=setdiff(vars,keep);
  if (~isscalar(vars)) CF=permute(CF,[i1,i2]); end;
  for i=length(vars):-1:length(s1)+1, CF=max(CF,[],i); end;
  msg=CF/sum(CF(:)); var=s1;    % set function outputs


