function dydt = dPlant2(t,y, Ff, Fh)
miu_max = 0.8384; KI = 24.3905; kd = 0.0209; kt = 0.1; kl = 0.7743;

Xv = y(1); Xd = y(2); Xl = y(3); Phi = y(4); V = 2; EFF = 0.0;

miuf = miu_max * 1 / ((Phi/KI).^3 + 1);
miud = kd + kt*Xl;


dXvdt = (miuf - miud - Ff/V + (1-EFF)*Fh/V)*Xv;
dXddt = miud*Xv - (kl + Ff/V - (1-EFF)*Fh/V)*Xd;
dXldt = kl*Xd - (Ff/V)*Xl;
dPhidt = Xv - (Ff/V)*Phi;

dydt = [dXvdt; dXddt; dXldt; dPhidt]; 
end