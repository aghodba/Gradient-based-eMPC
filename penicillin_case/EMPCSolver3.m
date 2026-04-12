function y = EMPCSolver3(xc,u,k,Concentrations,C,opt,NN,NPH,TS,sf,NCH,TF,Tint)


w=0.0;

options=optimoptions('fmincon','Algorithm','interior-point','Display','off','UseParallel','always',...
    'TolX',1e-6,'TolFun',1e-6,'TolCon',1e-6,'MaxFunctionEvaluations',50000,'MaxIterations',10000);
   
InitConc = [Concentrations.X(end) Concentrations.P(end) Concentrations.S(end)]; 

uEMPC = fmincon(@model_simulator_EMPC,ones(1,NCH)*u,[],[],[],[],ones(1,NCH)*0.2,ones(1,NCH)*12,[],options);
y = uEMPC;

    function Fp = model_simulator_EMPC(um)
        
        U = [ones(1,TF/Tint)*um(1) ones(1,TF/Tint)*um(2) ones(1,(NPH-2*NN)*TS/Tint + 1)*um(3)];
        input_FP = @(t) interp1(0:Tint:NPH*TS, U, t, 'linear', 'extrap');

        [t,MOutput] = ode15s(@(t,y) dModel(t,y,xc(1),xc(2),xc(3),xc(4),xc(5),[],xc(6),xc(7),[],sf,input_FP),(0:TS:NPH*TS),InitConc,opt);        
        
        Fp = - sum((MOutput(2:end,2) - C(end-NPH+1:end,2)).*input_FP(t(2:end))) + (w*(um-u).^2)*ones(NCH,1);        
    end
    
end    