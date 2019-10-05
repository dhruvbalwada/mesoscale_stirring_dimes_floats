function [ S2 ] = S2_from_E( E, k, r )

for i=1:length(r) 
    S2(i) = trapz(k, E.*(1- besselj(0, r(i)*k))); 
end

S2=(S2');

end

