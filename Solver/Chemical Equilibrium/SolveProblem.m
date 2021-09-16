function self = SolveProblem(self, i)
% Specify the problem type and the thermodynamic properties of the product
% mixture (P) required for the computations
% switch self.PD.ProblemType
%     case {'TP','TV'}
%         self.PS.strP{i} = SolveProblemTP_TV(self, self.PS.strR{i}, self.PD.pR.value, self.PD.TP.value);
%     case {'HP'}
%         if i==self.C.l_phi
%             self.PS.strP{i} = SolveProblemHP(self, self.PS.strR{i}, self.PD.phi.value(i), self.PD.pR.value);
%         else
%             self.PS.strP{i} = SolveProblemHP_EV_fast(self, self.PS.strR{i}, self.PD.pR.value,self.PS.strP{i+1});
%         end
%     case {'EV'}
%         if i==self.C.l_phi
%             self.PS.strP{i} = SolveProblemEV(self, self.PS.strR{i}, self.PD.phi.value(i), self.PD.pR.value);
%         else
%             self.PS.strP{i} = SolveProblemHP_EV_fast(self, self.PS.strR{i}, self.PD.pR.value,self.PS.strP{i+1});
%         end
%     case 'SP'
%         pP = self.PD.pP.value(i);
%         if i==self.C.l_phi
%             self.PS.strP{i} = SolveProblemSP(self, self.PS.strR{i}, self.PD.phi.value(i), pP);
%         else
%             self.PS.strP{i} = SolveProblemSP_fast(self, self.PS.strR{i}, self.PD.phi.value(i), pP, self.PS.strP{i+1});
%         end
%     case 'SV'
%         vP_vR = self.PD.vP_vR.value(i);
%         self.PS.strP{i} = SolveProblemSV(self, self.PS.strR{i}, self.PD.phi.value(i), vP_vR*self.PS.strR{i}.v);
%     case {'SHOCK_I','SHOCK_R'}
%         u1 = self.PD.u1.value(i);
%         if strcmp(self.PD.ProblemType,'SHOCK_I')
%             if i==self.C.l_phi
%                 [self.PS.strR{i}, self.PS.strP{i}] = SolveProblemSHOCK_I(self, self.PS.strR{i}, self.PD.phi.value(i), self.PD.pR.value, self.PD.TR.value, u1);
%             else
%                 [self.PS.strR{i}, self.PS.strP{i}] = SolveProblemSHOCK_I_fast(self, self.PS.strR{i}, self.PD.phi.value(i), self.PD.pR.value, self.PD.TR.value, u1, self.PS.strP{i+1});
%             end
%         else
%             if i==self.C.l_phi
%                 [self.PS.strR{i}, self.PS.str2{i}, self.PS.strP{i}] = SolveProblemSHOCK_R(self, self.PS.strR{i}, self.PD.phi.value(i), self.PD.pR.value, self.PD.TR.value, u1);
%             else
%                 [self.PS.strR{i}, self.PS.str2{i}, self.PS.strP{i}] = SolveProblemSHOCK_R_fast(self, self.PS.strR{i}, self.PD.phi.value(i), self.PD.pR.value, self.PD.TR.value, u1, self.PS.str2{i+1}, self.PS.strP{i+1});
%             end
%         end
%     case 'DET'
% %         self.PD.TR.value = self.PD.TR_vector.value(i);
% %         [self.PS.strR{i},self.PS.strP{i},self.TN.guess] = SolveProblemDET_main(self.PS.strR{i},self.PD.pR.value,self.PD.TR.value,self.TN.guess,self.E,self.S,self.C,self.M,self.PD,self.TN,self.strThProp);
%         [self.PS.strR{i},self.PS.strP{i},self.TN.guess] = SolveProblemDET_main_2(self, self.PS.strR{i},self.PD.pR.value,self.PD.TR.value, self.TN.guess);
% %         [self.PS.strR{i},self.PS.strP{i},self.TN.guess] = SolveProblemDET_hybrid_main(self.PS.strR{i},self.PD.pR.value,self.PD.TR.value,self.TN.guess,self.E,self.S,self.C,self.M,self.PD,self.TN,self.strThProp);
% %         [self.PS.strR{i},self.PS.strP{i},self.TN.guess] = SolveProblemDET_Jouget_main(self.PS.strR{i},self.PD.pR.value,self.PD.TR.value,self.TN.guess,self.E,self.S,self.C,self.M,self.PD,self.TN,self.strThProp);
%     case 'DET_OVERDRIVEN'
%         [self.PS.strR{i},self.PS.strP{i},self.TN.guess] = SolveProblemDET_OVERDRIVEN(self, self.PS.strR{i},self.PD.pR.value,self.PD.TR.value,self.TN.guess, self.PD.overdriven.value);
% end
% end

switch self.PD.ProblemType
    case {'SHOCK_I', 'SHOCK_R', 'DET', 'DET_OVERDRIVEN'}
        u1 = self.PD.u1.value(i);
        if strcmp(self.PD.ProblemType,'SHOCK_I')
            if i==self.C.l_phi
                [self.PS.strR{i}, self.PS.strP{i}] = SolveProblemSHOCK_I(self, self.PS.strR{i}, self.PD.phi.value(i), self.PD.pR.value, self.PD.TR.value, u1);
            else
                [self.PS.strR{i}, self.PS.strP{i}] = SolveProblemSHOCK_I_fast(self, self.PS.strR{i}, self.PD.phi.value(i), self.PD.pR.value, self.PD.TR.value, u1, self.PS.strP{i+1});
            end
        else
            if i==self.C.l_phi
                [self.PS.strR{i}, self.PS.str2{i}, self.PS.strP{i}] = SolveProblemSHOCK_R(self, self.PS.strR{i}, self.PD.phi.value(i), self.PD.pR.value, self.PD.TR.value, u1);
            else
                [self.PS.strR{i}, self.PS.str2{i}, self.PS.strP{i}] = SolveProblemSHOCK_R_fast(self, self.PS.strR{i}, self.PD.phi.value(i), self.PD.pR.value, self.PD.TR.value, u1, self.PS.str2{i+1}, self.PS.strP{i+1});
            end
        end
    otherwise
        if i == self.C.l_phi
            self.PS.strP{i} = equilibrate(self, self.PS.strR(i), self.PD.pP.value);
        else
            self.PS.strP{i} = equilibrate(self, self.PS.strR(i), self.PD.pP.value, self.PS.strP(i + 1));
        end
end
end