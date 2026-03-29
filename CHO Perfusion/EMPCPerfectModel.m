function OPT_T = EMPCPerfectModel
%LOAD DATA %opt = 3573
miu_max = 0.8384; KI = 24.3905; kd = 0.0209; kt = 0.1; kl = 0.7743;

err = 0.0;
err2 = 1 - 0.5*err;
NU=100;
NPH = 40;
Tint = 0.001; %time intervals in inputF
TS = 0.15; %Sampling time
TF = 0.3;
NN = TF/TS;
Ff0 = 1; %initial feedrate
Fh0 = 1;
NCH = 3;

InitConc = [0.2, 0, 0 ,0];

options=optimoptions('fmincon','Algorithm','interior-point','Display','off','UseParallel',false,...
    'TolX',1e-4,'TolFun',1e-4,'TolCon',1e-6,'MaxFunctionEvaluations',20000,'MaxIterations',5000);
options2 = optimoptions('patternsearch','MaxIterations',1000,'UseCompletePoll',true,'UseCompleteSearch',true,'UseParallel','always','PlotFcn','psplotbestf','FunctionTolerance',1e-4);
opt2=odeset('NonNegative',[1,2,3,4],'RelTol',1e-6,'AbsTol',1e-6, 'BDF', 'on');
opt=odeset('NonNegative',[1,2,3,4],'RelTol',1e-6,'AbsTol',1e-6, 'BDF', 'on');

u = ones(210,2).*[Ff0,Fh0];
UFf = ones(1,TF/Tint + 1)*Ff0;
UFh = ones(1,TF/Tint + 1)*Fh0;
input_Ff = @(t) interp1(0:Tint:TF, UFf, t, 'linear', 'extrap');
input_Fh = @(t) interp1(0:Tint:TF, UFh, t, 'linear', 'extrap');

Concentrations.Xv=[]; Concentrations.Xd=[]; Concentrations.Xl=[]; Concentrations.Phi=[]; Concentrations.XvNF=[];
[~,FirstPlantOutput] = ode15s(@(t,y) dPlant(t,y, miu_max, KI, kd, kt, kl, input_Ff, input_Fh), (0:TS:TF), InitConc, opt);
Concentrations.Xv = [FirstPlantOutput(1,1); FirstPlantOutput(2:3,1).*(err2 + err * rand(2,1))];
Concentrations.Xd = [FirstPlantOutput(1,2); FirstPlantOutput(2:3,2).*(err2 + err * rand(2,1))];
Concentrations.Xl = [FirstPlantOutput(1,3); FirstPlantOutput(2:3,3).*(err2 + err * rand(2,1))];
Concentrations.Phi = [FirstPlantOutput(1,4); FirstPlantOutput(2:3,4).*(err2 + err * rand(2,1))];
Concentrations.XvNF = [FirstPlantOutput(1,1); FirstPlantOutput(2:3,1)];

UTimes = [1];

%%

for k=2:NU
    
    u_EMPC = EMPCSolver_RealModel(u,k,Concentrations,opt,NPH,TS)

    u(k,:) = u_EMPC;
    
    UFf = [UFf ones(1,TF/Tint)*u(k,1)];
    UFh = [UFh ones(1,TF/Tint)*u(k,2)];
    input_Ff = @(t) interp1(0:Tint:TF*k, UFf, t, 'linear', 'extrap');
    input_Fh = @(t) interp1(0:Tint:TF*k, UFh, t, 'linear', 'extrap');
    
    [Xv1,Xd1,Xl1,Phi1] = plant_simulator();
    Concentrations.Xv =  [Concentrations.Xv; Xv1(NN*k:end).*(err2 + err * rand(NN,1))];
    Concentrations.Xd =  [Concentrations.Xd; Xd1(NN*k:end).*(err2 + err * rand(NN,1))];
    Concentrations.Xl =  [Concentrations.Xl; Xl1(NN*k:end).*(err2 + err * rand(NN,1))];
    Concentrations.Phi = [Concentrations.Phi; Phi1(NN*k:end).*(err2 + err * rand(NN,1))];
    Concentrations.XvNF = [Concentrations.XvNF; Xv1(NN*k:end)];
    
    
    display(k)
    save('Test5.mat')
    figure(2)
    subplot(2,2,1);plot(Concentrations.Xv);subplot(2,2,2);plot(u(1:k,1));subplot(2,2,3);plot(u(1:k,2));
    pause(0.01)
    
end


    function [Xvp,Xdp,Xlp,Phip] = plant_simulator()
        [~,POutput] = ode15s(@(t,y) dPlant(t,y, miu_max, KI, kd, kt, kl, input_Ff, input_Fh),(0:TS:k*TF),InitConc,opt2);
        Xvp = POutput(:,1);
        Xdp = POutput(:,2);
        Xlp = POutput(:,3);
        Phip = POutput(:,4);
    end

end