
variable x1, x2, yext, yabs, ysca;

require("ismdust");
(x1, x2) = linear_grid(1.0, 40.0, 10000);

fit_fun("ismdust(1, constant(1))");
yext = eval_fun(x1, x2);

fit_fun("ismdust_abs(1, constant(1))");
yabs = eval_fun(x1, x2);

fit_fun("ismdust_sca(1, constant(1))");
ysca = eval_fun(x1, x2);

% plots the cross-sections (tau, unitless)

xlog; ylog;
hplot(x1, x2, -log(yext), 1);
ohplot(x1, x2, -log(yabs), 2);
ohplot(x1, x2, -log(ysca), 4);
