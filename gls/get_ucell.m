
function uc = get_ucell(D, pN)

tD = D';

%[rn, un] = size(tDt);

uc = cell(pN,1);
[I,J,V] = find(tD);

for i = 1:pN
    uidx = (I == i);
    uc{i} = [J(uidx)'; V(uidx)'];
end






