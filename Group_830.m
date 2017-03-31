clear all
close all
clc
 
%% CVX exercise 
% To complete the exercise you need to fill in the missing part of the cvx
% proceedure (from cvx_begin to cvx_end) GOOD LUCK
 
L = 10; % Window size of the mpc proplem (Control horizon = Prediction horizon) %10
M = 200; % The duration of the control process %200
Ts = 1; % Time step (1 hour)
 
onesL = ones(L,1); % one vector of length L
zerosL = zeros(L,1); % zero vector of length L
 
% Defining vectors of the system parameters
E_A_sys = zeros(M-L+2,1);
Q_W_sys = zeros(M-L+1,1);
Q_G_sys = zeros(M-L+1,1);
Q_A_in_sys = zeros(M-L+1,1);
Q_A_out_sys = zeros(M-L+1,1);
Q_E_sys = zeros(M-L+1,1);
Q_bp_sys = zeros(M-L+1,1);
revenue = zeros(M-L+1,1);
 
% Generates the data (Low pass filtered white gaussian noise)
price_data_generator %(OBS: There is a rand seed = 1 in the file)
 
%% Plot price vectors
figure
hold on
stairs(P_E)
stairs(P_G,'r')
stairs(P_W,'g')
title('price data')
legend('Price Electricity','Price for gas','Price of burning waste')
ylabel('[DKK/MWh]')
xlabel('Sample [hour]')
%%
 
% Define the linear inequalities for the variables
Q_W_min = 0;
Q_W_max = 40;
Q_G_min = 0;
Q_G_max = 20;
E_A_min = 0;
E_A_max = 200;
Q_A_in_min = 0;
Q_A_in_max = 50;
Q_A_out_min = 0;
Q_A_out_max = 25;
 
E_A_sys(1) = 0; %Initial condition

for k = 1:M-L+1 % The main loop
 
cvx_begin % The begining of the optimization problem
 
% Define the variables %%% FILL IN %%%
    variables Q_E(L) Q_G(L) Q_W(L) E_A(L+1) Q_A_in(L) Q_A_out(L) Q_bp(L)  
    %E_A 
 
% Specify the optimization of cost %%% FILL IN %%% 
     %minimize( -(P_E(k:k+L-1)'*Q_E-(P_G(k:k+L-1)'*Q_G(1:L)+P_W(k:k+L-1)'*Q_W(1:L)))*Ts)
     maximize((P_E(k:L+k-1)'*Q_E-(P_G(k:L+k-1)'*Q_G + P_W(k:L+k-1)'*Q_W))*Ts)
    
% constraints %%% FILL IN %%%
    subject to 
       % Equations power flow:
       Q_W + Q_G == Q_bp + Q_A_in;
       Q_E == Q_bp + Q_A_out;
       % Save the first element of E_A which is calculated over a horizon
       E_A_sys(k) == E_A(1);
       % Constraints on the variables:
       Q_W_min*onesL <= Q_W <= Q_W_max*onesL; 
       Q_G_min*onesL <= Q_G <= Q_G_max*onesL; 
       E_A_min*ones(L+1,1) <= E_A <= E_A_max*ones(L+1,1);   
       Q_A_in_min*onesL <= Q_A_in <= Q_A_in_max*onesL; 
       Q_A_out_min*onesL <= Q_A_out <= Q_A_out_max*onesL; 
       % Accumulator Dynamics:
       E_A(2:L+1) == E_A(1:L)+ (Q_A_in - Q_A_out)*Ts;
       
 
cvx_end % The end of the optimization problem
cvx_status % Tells whether the problem is solved. 
%Does not tell you whether the problem is posed correctly. 

% Calculate the system (The first entry of the vectors)
% Save the data for analysis
E_A_sys(k+1) = E_A_sys(k) + (Q_A_in(1) - Q_A_out(1))*Ts; % The real system
Q_W_sys(k) = Q_W(1); 
Q_G_sys(k) = Q_G(1); 
Q_A_in_sys(k) = Q_A_in(1);
Q_A_out_sys(k) = Q_A_out(1);
Q_E_sys(k) = Q_E(1);
Q_bp_sys(k) = Q_bp(1);
revenue(k) = [-P_G(k); -P_W(k); P_E(k)]'*[Q_G(1); Q_W(1); Q_E(1)];
end
%% Plot the results
 
% I got you started! Make some more plots and investigate the results
close all
figure
stairs(E_A_sys)
title('The state of charge in the accumulator(E_Asys)')
ylabel('[MWh]')
xlabel('Sample [hour]')
 
figure
stairs(Q_W_sys)
title('Power from using waste(Q_Wsys)')
ylabel('[MW]')
xlabel('Sample [hour]')
 
figure
stairs(Q_G_sys)
title('Power from using gas(Q_Gsys)')
ylabel('[MW]')
xlabel('Sample [hour]')
 
figure
hold on
stairs(P_W,'r')
stairs(P_E)
stairs(P_G,'g')
legend('P_W','P_E','P_G')
title('Power prices')
ylabel('[DKK/MWh]')
xlabel('Sample [hour]')

figure
hold on
stairs(Q_A_in_sys,'r')
stairs(Q_A_out_sys)
legend('Q_A_{in}','Q_A_{out}')
title('Accumulator in- and outflow')
ylabel('[MW]')
xlabel('Sample [hour]')

figure
axis([0 120 -90 350])
hold on
grid on
stairs(P_G,'r','LineWidth',1)
stairs(P_E,'k','LineWidth',1)
stairs(Q_G_sys,'LineWidth',1)
legend('P_{Gas}','P_{Electricity}','Q_{Gas}')
title('Comparison of prices(P_G,P_W) with Q_G')
ylabel('[DKK & MW]')
xlabel('Sample [hour]')

figure
axis([0 120 -90 350])
hold on
grid on
stairs(P_W,'r','LineWidth',1)
stairs(P_E,'k','LineWidth',1)
stairs(Q_W_sys,'LineWidth',1)
legend('P_{Waste}','P_{Electricity}','Q_{Waste}')
title('Comparison of prices(P_E,P_W) with Q_W')
ylabel('[DKK & MW]')
xlabel('Sample [hour]')

figure
axis([0 120 -55 210])
hold on
grid on
stairs(E_A_sys,'r','LineWidth',1)
stairs(Q_E_sys,'LineWidth',1)
legend('E_A','Q_E')
title('Accumulator energy and output power')
ylabel('[J & MW]')
xlabel('Sample [hour]')