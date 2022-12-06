clear;
load mandrill;
colormap("gray");
[U,S,V] = svd(X);
image(U*S*V');

close all;
for i = 1:1:6 
    r = 2^i;

    dims = size(X);
    S_trunc = diag(S);
    S_trunc((r + 1):min(dims(1), dims(2))) = 0;
    S_trunc = diag(S_trunc);
    S_trunc(dims(1), dims(2)) = 0;

    gcf
    hold on;
    Xout = U * S_trunc * V';
    colormap("gray");
    image(flipud(Xout));
    title(strcat("r=", string(r)));
    saveas(gcf, strcat("DevamSisodraker_5a_", string(r), ".jpg"), "jpg");
    hold off;
end