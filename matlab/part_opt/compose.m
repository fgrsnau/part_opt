function P = compose(P1,P2)

A = to_class([1:size(P2,2)],class(P1));
P = mselect(P1,P2,repmat(A,size(P2,1),1));

end