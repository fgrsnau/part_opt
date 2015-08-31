function [V_R E_R B_R e_R] = restrict(G,R)

%% Lets rewrite thtis to compute outer boundary of R

ie = G.in(R,:); % this is sparse matrix too
oe = G.out(R,:); % this is sparse matrix

[ir ic iv] = find(ie);
[or oc ov] = find(oe);

V_R = unique([G.E(1,iv) G.E(2,ov)]);
B_R = setdiff(V_R,R);
[E_R m] = unique([G.E(:,iv) G.E(:,ov)]','rows');
E_R = E_R';
ee = [ic; oc]; 
e_R = ee(m);

n = length(V_R);
IR = zeros(1,G.nV);
IR(V_R) = 1:n;
E_R = IR(E_R);
B_R = IR(B_R);

return


n = length(R);
IR = zeros(1,G.nV);
IR(R) = 1:n;
% buld inverse R index, for R(1:n) find its number in R

ie = G.in(R,:); % this is sparse matrix too
oe = G.out(R,:); % this is sparse matrix

[ir ic iv] = find(ie);
[or oc ov] = find(oe);

% boundary: vertices in R with at least one in or out edge not in R
im = IR(G.E(1,iv))>0;
om = IR(G.E(2,ov))>0;
ibm = accumarray(ir,im>0,[n 1])<G.in_deg(R);
obm = accumarray(or,om>0,[n 1])<G.out_deg(R);
B_R = find(ibm | obm);

e_R = ov(om);
E_R = IR(G.E(:,ov(om)));% all edges R->R

%E_R = IR();
%mR(G.E(1,:))>0 & mR(G.E(2,:))>0;
%E_R = mR(G.E(:,mE_R));

end