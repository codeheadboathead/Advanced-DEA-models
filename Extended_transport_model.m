clc 
clear 
close
%Importing the excel files with data to matlab
O1 = readmatrix("output1.xlsx");
O2 = readmatrix("output2.csv.xlsx");
I = readmatrix("input .csv.xlsx");

%defining the two types of efficiencies
efficiency_bar_matrix = zeros(8,3);
efficiency_kink_matrix = zeros(8,3);


 %calulation of E bar matrix    
for i = 1:8
    for j = 1:3
     e = [-1*O1(i,j) -1*O2(i,j) 0];

     I_min_lhs = [0 0 I(i,j)];
     I_min_rhs = [1];
     
       A_inequal_lhs = [O1(i,1) O2(i,1) -1*I(i,1);O1(i,2) O2(i,2) -1*I(i,2);O1(i,3) O2(i,3) -1*I(i,3)];
       A_inequal_rhs = [0;0;0];
     
     lb = [0.000001 0.000001 0.000001];
     [x,eval] = linprog(e,A_inequal_lhs,A_inequal_rhs,I_min_lhs,I_min_rhs,lb) ;
     e_bar = -1*eval;
     efficiency_bar_matrix(i,j) = e_bar;
    end
end
    
%calculation for E kink martrix
   for i = 1:8
    for j = 1:3
     e = [-1*O1(i,j) -1*O2(i,j) 0];

     I_min_lhs = [0 0 I(i,j)];
     I_min_rhs = [1];
     
       A_inequal_lhs = [O1(1,j) O2(1,j) -1*I(1,j);O1(2,j) O2(2,j) -1*I(2,j);O1(3,j) O2(3,j) -1*I(3,j);
           O1(4,j) O2(4,j) -1*I(4,j);O1(5,j) O2(5,j) -1*I(5,j);
           O1(6,j) O2(6,j) -1*I(6,j);O1(7,j) O2(7,j) -1*I(7,j);O1(8,j) O2(8,j) -1*I(8,j);];
       A_inequal_rhs = [0;0;0;0;0;0;0;0];
     
     lb = [0.000001 0.000001 0.000001];
     [x,eval] = linprog(e,A_inequal_lhs,A_inequal_rhs,I_min_lhs,I_min_rhs,lb) ;
     e_kink = -1*eval;
     efficiency_kink_matrix(i,j) = e_kink;
    end
   end
 
   %calculation of E matrix from E bar and E kink
   efficiency = zeros(8,3);
   for i = 1:8
       for j = 1:3
           efficiency(i,j) = (efficiency_bar_matrix(i,j) + efficiency_kink_matrix(i,j))/2;
       end
   end

   %calculation of transportation plan
   ineff_matrix1x24 = zeros(1,24); 
   ineff_matrix8x3 = readmatrix("ineff.xlsx");

   %converiosn of 8x3 matrix into 1x24
   k=1;
   for i = 1:8
       for j = 1:3
           ineff_matrix1x24(k) = ineff_matrix8x3(i,j); 
           k = k+1;
       end
   end
   
   %importing the supply and demand matrices to matlab
supply = zeros(8,1);
c_to_s = readmatrix("complete_table.xlsx");

%creating the supply and demand vectors
for i = 1:8
   supply(i,1) = c_to_s(i,5);
end
demand = zeros(3,1);
for i = 1:3
   demand(i,1) = c_to_s(9,i+1);
end

%calculation of transportation plan from using inefficiency minimization

%Et = ineff_matrix8x3; %cost matrix
%s_equal_lhs = ones(8,3);
s_equal_rhs = supply; %supply vector
%d_equal_lhs = ones(3,8);
d_equal_rhs = demand'; %demand vector
%lbet = zeros(1,24);
% Convert the transportation problem into a linear programming problem
f = ineff_matrix1x24;  % Objective function coefficients
Aeq = readmatrix("Aeq_corrected.xlsx");  % Equality matrix
beq = [s_equal_rhs; d_equal_rhs(:)];  % Right-hand side of the equality matrix
lb = zeros(size(f))';  % Lower bounds of the decision variables

[x, fval] = linprog(f, [], [], Aeq, beq, lb, [], []);

% Reshape the solution into the transportation problem format
plan_matrix = zeros(8,3);
o=1;
for i = 1:8
    for j = 1:3
        plan_matrix(i,j) = x(o,1);
        o = o+1;

    end
end







    