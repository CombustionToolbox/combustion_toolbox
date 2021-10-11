function prop = read_CT(filename)
load(filename)
Nstruct = length(strP);
for j=Nstruct:-1:1
    prop.prop(1,j) = strP{1,j}.T;
    prop.prop(2,j) = strP{1,j}.rho;
    prop.prop(3,j) = strP{1,j}.S;
%     prop.prop(4,j) = strP{1,j}.G/strP{1,j}.mi;
    prop.prop(5,j) = strP{1,j}.h/strP{1,j}.mi;
    prop.prop(6,j) = strP{1,j}.e/strP{1,j}.mi;
    prop.prop(7,j) = strP{1,j}.p;
%     prop.prop(8,j) = strP{1,j}.u;
    prop.prop(9,j) = strP{1,j}.gamma;
    prop.prop(10,j) = strP{1,j}.hf/strP{1,j}.mi;
    prop.prop(11,j) = strP{1,j}.cP;
    prop.prop(12,j) = strP{1,j}.mi;
    prop.prop(13,j) = strP{1,j}.DhT/strP{1,j}.mi;
    prop.X(:,j) = strP{1,j}.Xi;
end
