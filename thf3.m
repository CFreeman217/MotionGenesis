[t_1,x_a] = ode23('doubleP',[0,100],[pi,pi,0,0]);
plot(t_1, x_a(:,3),t_1, x_a(:,4));
% ,t_1, x_a(:,2)
xlabel('t (s)');
ylabel('y (m)');
% plot(t_1, x_a(:,3),t_1, x_a(:,4));
% xlabel('t (s)');
% ylabel('y (m)')