Positive Variables A, B, PA, PB, XA, XB;
Variables Z;

Equations Objective, Unit_Price_A, Unit_Price_B, Production_A, Production_B, Limit;

Objective.. Z =E= A * PA + B * PB;
Unit_Price_A.. PA =E= 8 - A;
Unit_Price_B.. PB =E= 6 - 1.5*B;
Production_A.. A =E= 0.5 * XA;
Production_B.. B =E= 0.25 * XB;
Limit.. XA + XB =L= 10;

Model Max_Z /All/;
Solve Max_Z Using NLP Maximizing Z;

File Max_Z_Result /Max_Z_Result.txt/
Put Max_Z_Result
Put 'Results from Max Z' Put //
Put ' Water Available: 10 ', Put //
Put '          Units       Unit Price  Water Allocation ', Put /
Put 'A:', A.L, PA.L, XA.L, Put /
Put 'B:', B.L, PB.L, XB.L;
