SET t Time /1*15/;
PARAMETER Q(t) inflow
/
1        5
2        7
3        8
4        4
5        3
6        3
7        2
8        1
9        3
10       6
11       8
12       9
13       3
14       4
15       9
/;
SCALAR K /5/;
Variable Y;
POSITIVE VARIABLES S(t), R(t);

Equations Balance(t), Capacity(t);
Balance(t).. S(t+1)$(ORD(t) lt 15) + S('1')$(ORD(t) eq 15) =E= S(t) + Q(t) - Y - R(t);
Capacity(t).. S(t) =L= K;

Model Max_Y / ALL /;
Solve Max_Y Using LP Maximizing Y;

File Max_Y_Result /Max_Y_Result.txt/
Put Max_Y_Result
Put 'Results from Max Y' Put //
Put ' Yield ', Put Y.L, Put /
Put ' Capacity ', Put K, Put //
Put ' t                  Inflow      Storage     Release ' Put /
Put '                    Q(t)        S(t)        R(t)    ' Put /
LOOP(t, Put t.TL, Q(t), S.L(t), R.L(t), Y.L/)
