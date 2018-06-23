
variable x1, x2, y;

require("olivineabs");
(x1, x2) = linear_grid(1.0, 2.5, 1000);
fit_fun("olivineabs(1, powerlaw(1))");
set_par(3, 100.0);
y = eval_fun(x1, x2);
xlog; ylog;
hplot(x1, x2, y);
