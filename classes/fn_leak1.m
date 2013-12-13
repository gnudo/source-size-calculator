function n=fn_leak1(x,p)
% leakage function by R. Mokso
% x is the array, p is the approximate period

for k=1:1:round(p)
    t=x(1:end-k+1);
    ft=fft(t);
    aft=abs(ft);
    naft=aft/(length(t));
    pos=round(length(t)/p);
    [peak,peakpos]=max(naft(pos-2:pos+2));
    peakpos=pos-2+peakpos-1;
    %disp([peak*1e3])
    f1(k)=peak;
end

%plot(f1);

[peak,peakpos]=max(f1);
n=length(x)-peakpos+1;