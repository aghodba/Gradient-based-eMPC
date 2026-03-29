function y = EMPCSolver_RealModel(u,k,Concentrations,opt,NPH,TS)
NU=100;
TF = 0.3;
options=optimoptions('fmincon','Algorithm','interior-point','Display','off','UseParallel','always',...
    'TolX',1e-6,'TolFun',1e-6,'TolCon',1e-6,'MaxFunctionEvaluations',50000,'MaxIterations',10000);

InitConc = [Concentrations.Xv(end) Concentrations.Xd(end) Concentrations.Xl(end) Concentrations.Phi(end)];


% [uEMPC,~,ef] = fmincon(@model_simulator_EMPC,[1,0.9],[],[],[],[],[0 0],[2 2],@nlin,options)

ms = MultiStart;
problem = createOptimProblem('fmincon','x0',[1,0.9],'objective',@model_simulator_EMPC,'lb',[0,0],'ub',[10,10], 'nonlcon',@nlin);
uEMPC = run(ms,problem,10)

y = uEMPC;

    function Fp = model_simulator_EMPC(um)
        
        [~,MOutput] = ode15s(@(t,y) dPlant2(t,y, um(1), um(2)),(0:TS:(NU-k+1)*TF),InitConc,opt);
        
        Fp = - (sum(MOutput(2:end,1))*um(2) + MOutput(end,1)*2);
        
    end

    function [c,ceq] = nlin(um)
        [~,MOutput1] = ode15s(@(t,y) dPlant2(t,y, um(1), um(2)),(0:TS:(NU-k+1)*TF),InitConc,opt);
        c = [um(2)-um(1)
            0.7 - (MOutput1(end,1)./(MOutput1(end,1)+MOutput1(end,2)))];
        ceq = [];
    end

end