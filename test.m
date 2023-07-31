% -----------------
% test 
l_D = l_bound_D(1);
u_D = u_bound_D(1);
l_L = l_bound_L(1);
u_L = u_bound_L(1);

x_opt = [((u_D-l_D).*rand(3,1)+l_D)',((u_L-l_L).*rand(3,1)+l_L)'];
% -----------------