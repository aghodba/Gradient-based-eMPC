function dydt = CHOModel(t,y)
miu_max = 0.8384; KI = 24.3905; kd = 0.0209; kt = 0.0290; kl = 0.7743;

Xv = y(1); Xd = y(2); Xl = y(3); Phi = y(4); V = y(5);

Ff=10;
Fh=9.99;
Fb=0.01;

miuf = miu_max * 1 / ((Phi/KI).^3 + 1);
miud = kd + kt*Xl;

dXvdt = (miuf - miud - Ff/V + Fh/V)*Xv;
dXddt = miud*Xv - (kl + Ff/V - Fh/V)*Xd;
dXldt = kl*Xd - (Ff/V)*Xl;
dPhidt = Xv - (Ff/V)*Phi;
dVdt = Ff - Fh - Fb;

dydt = [dXvdt; dXddt; dXldt; dPhidt; dVdt]; 
end