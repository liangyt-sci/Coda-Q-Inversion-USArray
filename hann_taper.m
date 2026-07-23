function arrout = hann_taper(arrin, delta, npts, len_seconds)
%HANN_TAPER Apply a half-Hann taper to both ends of a time series.

nlen = floor(len_seconds / delta);
arrout = arrin(:);

if nlen <= 1 || 2*nlen > npts
    return;
end

arg = (0:nlen-1) * pi / (2*nlen);
taper = sin(arg(:));

arrout(1:nlen) = arrout(1:nlen) .* taper;
arrout(npts-nlen+1:npts) = arrout(npts-nlen+1:npts) .* flipud(taper);
end
