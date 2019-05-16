function [fcn] = doubleP(x, t)
    fcn = zeros(4,1);
    l1 = 1;
    l2 = 2;
    m1 = 2;
    m2 = 2;
    g = 9.8;
    u_1 = t(1);
    v_1 = t(2);
    u_2 = t(3);
    v_2 = t(4);
    a = (m1 + m2) * l1;
    b = m2*l2*cos(u_1-v_1);
    c = - m2*l2*v_2^2*sin(u_1-v_1) - g*(m1+m2)*sin(u_1);
    d = m2*l2;
    e = m2*l1*cos(u_1-v_1);
    f = m2*l1*u_2^2*sin(u_1-v_1) - m2*g*sin(v_1);
    fcn(1) = u_1;
    fcn(2) = v_1;
    fcn(3) = (c*d - b*f) / (d*e - d*a);
%     fcn(4) = (d*e - f*b) / (b*c - d*a);
    fcn(4) = ((a*b*f - c*e) / (b*(d*a - e)));
end