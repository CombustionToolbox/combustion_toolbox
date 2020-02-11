function prop = read_CT_strR(filename)
load(filename)
Nstruct = length(strR);
for j=Nstruct:-1:1
    prop.prop(1,j) = strR{1,j}.T;
    prop.prop(2,j) = strR{1,j}.rho;
    prop.prop(3,j) = strR{1,j}.S;
%     prop.prop(4,j) = strR{1,j}.G/strR{1,j}.mi;
    prop.prop(5,j) = strR{1,j}.h/strR{1,j}.mi;
    prop.prop(6,j) = strR{1,j}.e/strR{1,j}.mi;
    prop.prop(7,j) = strR{1,j}.p;
    prop.prop(8,j) = strR{1,j}.u;
    prop.prop(9,j) = strR{1,j}.gamma;
    prop.prop(10,j) = strR{1,j}.hf/strR{1,j}.mi;
    prop.prop(11,j) = strR{1,j}.cP;
    prop.prop(12,j) = strR{1,j}.mi;
    prop.prop(13,j) = strR{1,j}.DhT/strR{1,j}.mi;
    prop.X(:,j) = strR{1,j}.Xi;
end
