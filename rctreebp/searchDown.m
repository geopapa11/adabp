function T=searchDown(T,v)
  vars = updateMAP(T,v);
  if (any(T.diff(vars)))			% found a change?
    %fprintf('searchDown change at %d\n',v);
    T=searchChange(T,v);			%  check for cascading changes
  end;
  return;

function T=searchChange(T,v)			% called when a change is discovered at v
  children=T.Id(T.Scars{T.Latest(v)}); 		% Get children of v in the tree
  for u=children				% check each of them for potential changes
    T=checkForChanges(T,u,v);
  end;

% 
function T=checkForChanges(T,u,v)
  %fprintf('Checking %d from %d : ',u,v);
  vars = updateMAP(T,u);			% compute and update MAP for u

  if (any(T.diff(vars)))			% if any of them have changed
    %fprintf('changed\n');
    T=searchChange(T,u);			%  we need to search u's children as well
  else						% if not, continue checking down and toward v:
    %fprintf('nope\n');
    nextchildren = T.Id(T.Scars{T.Latest(u)});	%  get u's children

    for w=nextchildren,				%  find any that were between u and v & check them.
      %orient(T,w,0);				% re-orient w just in case (?)
      if (T.in(w)==v || T.out(w)==v) T=checkForChanges(T,w,v); end;
    end;
  end;

