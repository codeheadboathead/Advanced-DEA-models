clc
clear
close
det_input = [5 8 7];
stochastic_input_mu = [14 15 12];
stochastic_input_var = [1.4 1.5 1.2];
covariance_martix = [0 0.9 0.6; 0.9 0 0.7; 0.6 0.7 0];
output_diploma = [9 5 4];
output_bacholers = [4 7 9];
output_masters = [16 10 13];
complete_output = [output_diploma; output_bacholers; output_masters];
lb = [0 0 0 0];

obejctive_function = [1 0 0 0];%[phi0 lamda1 lamda2 lamda3]

diff_matrix = zeros(1,3);
add_matrix = zeros(1,3);

for i = 1:3
    diff_matrix(1,i) = stochastic_input_mu(1,i) - 1.96*sqrt(stochastic_input_var(1,i));
end
for i = 1:3
    add_matrix(1,i) = stochastic_input_mu(1,i) + 1.96*sqrt(stochastic_input_var(1,i));
end

efficecny_matrix = zeros(1,3);

for i = 1:3
linineq_stochastic_input_lhs =  [-1*diff_matrix(1,i) diff_matrix(1,1) diff_matrix(1,2) diff_matrix(1,3);-1*add_matrix(1,i) add_matrix(1,1) add_matrix(1,2) add_matrix(1,3);-det_input(1,i) det_input;0 -1*output_diploma;0 -1*output_bacholers;0 -1*output_masters];
linineq_stochastic_input_rhs = [0 ; 0 ; 0 ; -output_diploma(1,i); -output_bacholers(1,i);-output_masters(1,i)];

lineq_BCC_cond_lhs = [0 1 1 1];
lineq_BCC_cond_rhs = [1];
lb = [0 0 0 0];

[f,fval] = linprog(obejctive_function,linineq_stochastic_input_lhs,linineq_stochastic_input_rhs,lineq_BCC_cond_lhs,lineq_BCC_cond_rhs,lb);

efficecny_matrix(1,i) = fval;
end






