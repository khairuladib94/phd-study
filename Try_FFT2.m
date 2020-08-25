 t = 0:0.001:1-0.001;
 x = cos(2*pi*100*t)+cos(2*pi*200*t);
 plot(x)
 xdft = fft(x);
 % the DFT bins for 100 Hz are 101 and 1000-101+2
 indices100 = [101 length(xdft)-101+2];
 ydft = zeros(size(xdft));
 ydft(indices100) = xdft(indices100);
 y = ifft(ydft,'symmetric');
 
 figure
 plot(y)
