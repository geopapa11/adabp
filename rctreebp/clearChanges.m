function T=clearChanges(T)
  for i=1:T.changedCount,
   modU32Array(T.diff,T.changed(i),uint32(0));
   modU32Array(T.changed,i,uint32(0));
  end;
  modU32Array(T.changedCount,1,uint32(0));


