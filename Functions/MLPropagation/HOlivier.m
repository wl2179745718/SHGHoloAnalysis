function H = HOlivier(Vn,X,Us,AG)
%The H function defined by Olivier on page 10 of his note

H = myifft((myfft((Us + X).*Vn)).*(AG));

end