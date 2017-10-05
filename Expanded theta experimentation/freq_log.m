function [freq]=freq_log(freq_0,freq_end,number)
step=(log(freq_end)-log(freq_0))/(number-1);
for n=1:number
    f=exp(log(freq_0)+step*(n-1));
    freq(n)=f;
end
