function dydt = dPlant(t,y,mx,Kx,mp,Kp,KI,KH,Yxs,Yps,Mx,sf,input_FP)

F = input_FP(t);
dydt(1,1) = mx*y(3)*y(1)/(Kx*y(1)+y(3)) - F*y(1)/120;
dydt(2,1) = mp*y(3)*y(1)/(Kp+y(3)+(y(3)^2)/KI) - KH*y(2) - F*y(2)/120;
dydt(3,1) = (-1/Yxs)*(mx*y(3)*y(1)/(Kx*y(1)+y(3))) - (1/Yps)*(mp*y(3)*y(1)/(Kp+y(3)+(y(3)^2)/KI)) - Mx*y(1) + F*sf/120 - F*y(3)/120;
end