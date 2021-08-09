function phi_c = check_phi_c(strR,pP,TP,E,S,C,M,PD,TN,strThProp)
    Fuel = PD.Fuel;
    A0 = C.A0.Value;
    [~,~,~,~,~,~,~,~,~,phi_c,~] = CalculateProductsCC(strR.NatomE,1,TP,E.Elements,TN.factor_c,PD.Fuel,strThProp);
    phi = 0.99*phi_c;
    PD.S_Oxidizer = {'O2'}; PD.N_Oxidizer = PD.phi_t/phi;
    PD.R_Oxidizer = SetSpecies(C.M0.Value, PD.S_Oxidizer, PD.N_Oxidizer, PD.TR.Value,find_idx( PD.S_Oxidizer, S.NameSpecies),strThProp);
    PD.S_Inert = {'N2'};  PD.N_Inert =  PD.phi_t/phi*79/21;
    PD.R_Inert = SetSpecies( C.M0.Value, PD.S_Inert, PD.N_Inert, PD.TR.Value,find_idx( PD.S_Inert, S.NameSpecies),strThProp);
    R =  PD.R_Fuel+ PD.R_Oxidizer+ PD.R_Inert;
    strR = ComputeProperties(C.A0.Value,R, PD.pR.Value, PD.TR.Value, E.ind_C, E.ind_H);

    [NCO2P_0,NCOP_0,NH2OP_0,NH2P_0,NO2P_0,NN2P_0,NHeP_0,NArP_0,NCgrP_0,phi_c,FLAG_SOOT] = CalculateProductsCC(strR.NatomE,phi,TP,E.Elements,TN.factor_c,PD.Fuel,strThProp);
    P = SetSpecies(C.M0.Value,S.List_Compute_Species,[NCO2P_0,NCOP_0,NH2OP_0,NH2P_0,NO2P_0,NN2P_0,NHeP_0,NArP_0,NCgrP_0]',TP,S.idx_fixed,strThProp);
    [P, ~] = CalculateProductsIC(P,TP,pP,strR.v,phi,M.minor_products,phi_c,FLAG_SOOT,C.M0.Value,C.A0.Value,E,S,C,M,PD,TN,strThProp);
    if P(S.idx_CO2,1) >= 1e-2 || P(S.idx_H2O,1) >= 1e-2
        phi_c = (2*(Fuel.x+Fuel.y/4-Fuel.z/2))/(sum(P(:,1).*A0(:,E.ind_O)) - Fuel.z); % C_x H_y O_z
    end
end