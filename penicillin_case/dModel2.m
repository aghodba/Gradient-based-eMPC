function dydt = dModel2(t,y,flag,K,input_FP,extra_var)
dydt=zeros(3,1);

mx=K(1);
Kx=K(2);
mp=K(3);
Kp=K(4);
KI=K(5);
Yxs=K(6);
Yps=K(7);

sf = extra_var(1);


F = input_FP(t);
dydt(1,1) = mx*y(3)*y(1)/(Kx*y(1)+y(3)) - F*y(1)/120;
dydt(2,1) = mp*y(3)*y(1)/(Kp+y(3)+(y(3)^2)/KI) - F*y(2)/120;
dydt(3,1) = (-1/Yxs)*(mx*y(3)*y(1)/(Kx*y(1)+y(3))) - (1/Yps)*(mp*y(3)*y(1)/(Kp+y(3)+(y(3)^2)/KI)) + F*sf/120 - F*y(3)/120;

end