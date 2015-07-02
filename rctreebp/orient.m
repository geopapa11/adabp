function orient(T,v,recurse)
  if (T.in(v)==0) return; end;
  u=T.in(v); w=T.out(v);
  if (recurse)
    if (w==0 || T.in(u)==w || T.out(u)==w) orient(T,u,recurse); else orient(T,w,recurse); end;
  end;
  %re-orient v given orientations of parents:
  if (w~=0 && (T.in(u)==w || T.out(w)==u)),
    modU32Array(T.in,v,w); modU32Array(T.out,v,u);  %T.in(v)=w; T.out(v)=u;
  else
    modU32Array(T.in,v,u); modU32Array(T.out,v,w);  %T.in(v)=u; T.out(v)=w;
  end;


  
