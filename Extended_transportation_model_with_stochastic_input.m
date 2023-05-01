%Extended transportation model with stochastic input 
clc 
close
clear

O1 = readmatrix("output1.xlsx");
O2 = readmatrix("output2.csv.xlsx");
I_raw = readmatrix("input_data_stochastic.xlsx");

I_stochastic_mean = readmatrix("input .csv.xlsx");
I_stochastic_var = zeros(8,3);
k = 1;
for i = 1:8
    for j = 1:3
    
    I_stochastic_var(i,j) = I_raw(k,2);
    k = k+1;
    end
end

diff_matrix = zeros(1,3);
add_matrix = zeros(1,3);

for i = 1:8
    for j = 1:3
    diff_matrix(i,j) = I_stochastic_mean(i,j) - 1.96*sqrt(I_stochastic_var(i,j));
    add_matrix(i,j) = I_stochastic_mean(i,j) + 1.96*sqrt(I_stochastic_var(i,j));
    
    end
end


efficiency_bar_matrix = zeros(8,3);
efficiency_kink_matrix = zeros(8,3);

 %calulation of E bar matrix    
 for i=1:8
     for j = 1:3
 e_bar = [1 0 0 0];
 linineq_stochastic_input_lhs_bar =  [-1*diff_matrix(i,j) diff_matrix(i,1) diff_matrix(i,2) diff_matrix(i,3);-1*add_matrix(i,j) add_matrix(i,1) add_matrix(i,2) add_matrix(i,3);0 -1*O1(i,:);0 -1*O2(i,:)];
 linineq_stochastic_input_rhs_bar = [0 ; 0 ; -O1(i,j); -O2(i,j)];

 lb = [0 0 0 0];
 lineq_BCC_cond_lhs = [0 1 1 1];
 lineq_BCC_cond_rhs = [1];

[fbar,fvalbar] = linprog(e_bar,linineq_stochastic_input_lhs_bar,linineq_stochastic_input_rhs_bar,lineq_BCC_cond_lhs,lineq_BCC_cond_rhs,lb);

efficiency_bar_matrix(i,j) = fvalbar;
     end
 end

 %calulation of E kink matrix    
 for i=1:8
     for j = 1:3
 e_kink = [1 0 0 0 0 0 0 0 0];
 linineq_stochastic_input_lhs_kink =  [-1*diff_matrix(i,j) diff_matrix(1,j) diff_matrix(2,j) diff_matrix(3,j) diff_matrix(4,j) diff_matrix(5,j) diff_matrix(6,j) diff_matrix(7,j) diff_matrix(8,j);-1*add_matrix(i,j) add_matrix(1,j) add_matrix(2,j) add_matrix(3,j) add_matrix(4,j) add_matrix(5,j) add_matrix(6,j) add_matrix(7,j) add_matrix(8,j);0 -1*O1(:,j)';0 -1*O2(:,j)'];
 linineq_stochastic_input_rhs_kink = [0 ; 0 ; -O1(i,j); -O2(i,j)];

 lb = [0 0 0 0 0 0 0 0 0];
 lineq_BCC_cond_lhs = [0 1 1 1 1 1 1 1 1];
 lineq_BCC_cond_rhs = [1];

[fkink,fvalkink] = linprog(e_kink,linineq_stochastic_input_lhs_kink,linineq_stochastic_input_rhs_kink,lineq_BCC_cond_lhs,lineq_BCC_cond_rhs,lb);

efficiency_kink_matrix(i,j) = fvalkink;
     end
 end

  %calculation of E matrix from E bar and E kink
   efficiency = zeros(8,3);
   for i = 1:8
       for j = 1:3
           efficiency(i,j) = (efficiency_bar_matrix(i,j) + efficiency_kink_matrix(i,j))/2;
       end
   end

   %calculation of ineff matrix from efficinecy
   ineff_matrix8x3 = zeros(8,3);

   for i=1:8
       for j =1:3
       ineff_matrix8x3(i,j) = 1 - efficiency(i,j);
       end
   end

      ineff_matrix1x24 = zeros(1,24);

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


