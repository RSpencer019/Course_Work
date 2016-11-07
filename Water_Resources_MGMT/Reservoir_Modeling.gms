SETS t / Jan, Feb, Mar, Apr, May, Jun, Jul, Aug, Sep, Oct, Nov, Dec /;
SET crop /Cotton, Wheat, Rice, Lucerne/;

Parameters Q(t)
/
Jan      386
Feb      346
Mar      416
Apr      713
May      1532
Jun      2373
Jul      2157
Aug      1431
Sep      798
Oct      582
Nov      487
Dec      420
/

G(crop)
/
Cotton   350
Wheat    70
Rice     225
Lucerne  50
/

A1Max(crop)
/
Cotton   500000
Wheat    500000
Rice     0
Lucerne  50000
/

A2Max(crop)
/
Cotton   65000
Wheat    30000
Rice     500000
Lucerne  500000
/;

TABLE water_dem(crop,t) crop water demands (mm per month)
         Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec
Cotton   150 150 100 100 200 250 250 200 50  0   0   0
Wheat    0   0   0   100 300 300 0   0   200 0   0   0
Rice     0   0   0   450 700 800 800 700 50  0   0   0
Lucerne  0   0   0   50  200 350 350 250 50  0   0   0
;

Variable Z;

Positive Variables W1(crop,t), W2(crop,t), D1(t), D2(t), D1(t), D2(t), R1(t), R2(t), R3(t), R4(t), K(t), A1(crop), A2(crop);

Equations Area_Country1, Area_Country2, WaterDemand_Country1(crop,t), WaterDemand_Country2(crop,t), Objective(crop,t), R1eqn(t), R2eqn(t), R3eqn(t), R4eqn(t), Keqn1, Keqn2, Keqn3, Keqn4;

Keqn1.. K('Jan') =E= 14000;
Keqn2(t).. K(t) =G= 5500;
Keqn3(t).. K('Dec') =E= 13000;
Keqn4(t).. K(t+1)$(ORD(t) ne 11) + K('Jan')$(ORD(t) eq 11) =E= K(t) + Q(t) - R1(t);

R1eqn(t).. R1(t) =L= K(t)+Q(t)-5500;
R2eqn(t).. R2(t) =E= R1(t) - D1(t);
R3eqn(t).. R3(t) =E= R2(t) + (0.5*D1(t)) - D2(t);
R4eqn(t).. R4(t) =E= R3(t) + (0.5*D2(t));

Area_Country1(crop).. A1(crop) =L= A1Max(crop);
Area_Country2(crop).. A2(crop) =L= A2Max(crop);

WaterDemand_Country1(crop,t).. W1(crop, t)*A1(crop)/100000 =L= D1(t);
WaterDemand_Country2(crop,t).. W2(crop, t)*A2(crop)/100000 =L= D2(t);

alias (crop, ai);
Objective(crop, t).. Z =E= SUM(ai,G(ai)*(A1(ai)+A2(ai)));

Model Max_Z /All/;
Solve Max_Z Using NLP Maximizing Z;
