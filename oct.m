function out = oct(CR,SR)

out = 1;
if CR < 0
    out = out + 4;
end

if SR < 0
    out = out + 2;
end

if abs(CR) - abs(SR) < 0
    out = out + 1;
end

