function [A_normalized] = normalize(A, B, thrsh)
%UNTITLED3 Summary of this function goes here

cond = abs(A) > thrsh & abs(B) > thrsh;

A_select = A.*cond;
B_select = B.*cond;

fun = @(x)max(max(abs((x(1)+1i*x(2))*A_select - B_select)));
x0 = [1,0];
x_factor = fminsearch(fun,x0);

A_normalized = (x_factor(1)+1i*x_factor(2))*A;

end