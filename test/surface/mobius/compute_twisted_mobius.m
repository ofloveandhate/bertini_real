% [fn] = compute_twisted_mobius(twist)
%
% fn [out] : a symbolic function of x, y, z and parameters a, b
% twist [in] : an integer describing the twist of the surface
%
% Created by Travis C. Wert at the University of Notre Dame under Dani Brake
% on 4.26.17
%
% this work is inspired by this paper:
% https://imaginary.org/sites/default/files/moebiusband.pdf

function fn = compute_twisted_mobius(twist)



syms x y z t a b c2 s2 two_cs c s C S left1 left1a left1b left2 right1 right1a right1b right2

E_psi = c2*(a*(t-1)^2+b*z^2)+two_cs*(a-b)*(t-1)*z+s2*(b*(t-1)^2+a*z^2)-a*b; %from original paper found here https://imaginary.org/sites/default/files/moebiusband.pdf

assume(C,'real'); %assists with the substitutions later
assume(S,'real');
left1 = real(expand((C+1i*S)^twist)); %derived from De Moivre's formula
left2 = imag(expand((C+1i*S)^twist));

%There will be three equations that we need, c^2, s^2, and 2cs
%right1a will get s^2, right1b will get c^2, and right2 will get 2cs

right1 = c^2 - s^2;
right2 = 2*c*s;

right1a = subs(right1,c^2,1-s^2);
left1a = -(left1-1)/2;
left1a = subs(left1a,C,x/t);
left1a = subs(left1a,S,y/t);
left1a = expand(left1a);

right1b = subs(right1,s^2,1-c^2);
left1b = (left1+1)/2;
left1b = subs(left1b,C,x/t);
left1b = subs(left1b,S,y/t);
left1b = expand(left1b);

left2 = subs(left2,S,y/t);
left2 = subs(left2,C,x/t);

%substitute each of the left hand sides into E_psi
E_psi = subs(E_psi,{c2,two_cs,s2}, ...
    {left1b, ...
    left2, ...
    left1a});

E_psi = (E_psi)* (2*t^twist); %cancel the demonimator for the equation
E_psi = expand(E_psi);
simplify(E_psi);

%gather the odd and even coeffecients of t
[ct,tt] = coeffs(E_psi,t);

even = sym(0); %seed the loop
odd = sym(0);

for ii = 1:length(tt)
   deg = feval(symengine,'polylib:degree',char(tt(ii)),'t');
   if mod(deg,2) == 1 %odd
       odd = odd + ct(ii)*tt(ii);
   else % even
       even = even + ct(ii)*tt(ii);
   end
end

assume(a>0 | a<1);
assume(b>0 | b<1);
odd = collect(odd/t,t); %took out a t term, to be added back later
even = collect(even,t);
even = simplify(even);

%square both sides
even = even*even;
odd = odd*odd;

even = subs(even,t,sqrt(x^2+y^2));
odd = subs(odd,t,sqrt(x^2+y^2))*(x^2+y^2); %here the original t term (now t^2) is added back

fn = even-odd;
end
