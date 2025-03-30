\ENCODING=ISO-8859-1
\Problem name: ilog.cplex

Maximize
 obj1: 2 p_0 + p_1 + 1.5 p_2 + 2.5 p_3 + 0.5 p_4
Subject To
 c1:  p_4 - c_2 >= 0
 c2:  4 p_0 + 2 p_1 + 3 p_2 + 5 p_3 + p_4 >= 11
 c3:  4 p_0 + 2 p_1 + 3 p_2 + 5 p_3 + p_4 >= 5
 c4:  4 p_0 + 2 p_1 + 3 p_2 + 5 p_3 + p_4 <= 12
 c5:  3 p_0 + p_3 - 2 c_0 - 2 c_1 - 2 c_3 <= 0
 c6:  p_1 + p_4 - 2 c_2 - c_0 - c_1 - c_3 - c_4 <= 0
 c7:  p_0 + p_2 + 2 p_3 - c_0 - 2 c_1 - 2 c_4 <= 0
 c8:  p_1 + p_3 - c_2 - c_3 - c_4 <= 0
 c9:  2 p_2 + p_3 - 2 c_2 - c_0 - c_1 - c_3 - 2 c_4 <= 0
 c10: c_2 + c_0 + c_1 + c_3 + c_4  = 2
Bounds
 0 <= p_0 <= 1
 0 <= p_1 <= 1
 0 <= p_2 <= 1
 0 <= p_3 <= 1
 0 <= p_4 <= 1
 0 <= c_2 <= 1
 0 <= c_0 <= 1
 0 <= c_1 <= 1
 0 <= c_3 <= 1
 0 <= c_4 <= 1
Binaries
 p_0  p_1  p_2  p_3  p_4  c_2  c_0  c_1  c_3  c_4 
End
