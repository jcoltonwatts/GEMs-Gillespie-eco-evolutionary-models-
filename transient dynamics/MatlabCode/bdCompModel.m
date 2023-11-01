function dydt = bdCompModel(~,y,b_max1,b11, b_max2, b22, d_min1, d11, d_min2, d22, b12, b21, d12, d21)
dydt = zeros(size(y));

% variables
R1 = y(1);
R2 = y(2);

dydt(1) = (b_max1 - b11*R1 - b12*R2)*R1 - (d_min1 + d11*R1 + d12*R2)*R1;
dydt(2) = (b_max2 - b22*R2 - b21*R1)*R2 - (d_min2 + d22*R2 + d21*R1)*R2;