% Triangle pulse implementation
function res = triangle_pulse(t,w,h)
    t = t*(2/w);
    res = 1-min(1*abs(t-1),1);
    res = h*res;
end 