% Triangle pulse derivative 
function res = triangle_pulse_derivative(t,w,h)
    t = t*(2/w);
    res = -sign(t-1).*(abs(t-1)<1);
    res = (2*h/w)*res;
end