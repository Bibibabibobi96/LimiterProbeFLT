function [B_next] = interp3_nearest_Liao(R, Z, phi, B, r, z, t,R_min,R_max,Z_min,Z_max)

% global n nz nr nphi
if sum(isnan([r, z, t])) == 0

    nr=(r-R_min)/((R_max-R_min)/(length(R)-1))+1;
    nr_pre = floor(nr);
    nr_next = nr_pre+1;

    nz=(z-Z_min)/((Z_max-Z_min)/(length(Z)-1))+1;
    nz_pre = floor(nz);
    nz_next = nz_pre+1;


    nt=(t)/((2*pi)/(length(phi)-1))+1;
    nt_pre = floor(nt);
    nt_next = nt_pre+1;



    while  nt_pre > length(phi) + 1

        nt_pre = nt_pre-length(phi);
        nt_next = nt_next-length(phi);
    end
    if nt_pre == length(phi)+1

        nt_next = nt_next-length(phi);

    end

    while  nt_pre<0

        nt_pre = nt_pre+length(phi);
        nt_next = nt_next+length(phi);
    end

    if   nt_pre==0
        nt_pre = length(phi);

    end



    rd=nr-nr_pre;
    zd=nz-nz_pre;
    td=nt-nt_pre;
    % for i=1:length(r)
    %     if nr_next(i)>=350
    %      [B000(i),B100(i),B010(i),B001(i),B101(i),B011(i),B110(i),B111(i)] = deal(NaN);
    %     elseif nz_next(i)>=600
    % %      B000(i) = NaN;
    %      [B000(i),B100(i),B010(i),B001(i),B101(i),B011(i),B110(i),B111(i)] = deal(NaN);
    %     else
    %      B000(i) = B(nr_pre(i),nz_pre(i),nt_pre(i));
    %      B100(i) = B(nr_next(i),nz_pre(i),nt_pre(i));
    %      B010(i) = B(nr_pre(i),nz_next(i),nt_pre(i));
    %      B001(i) = B(nr_pre(i),nz_pre(i),nt_next(i));
    %      B101(i) = B(nr_next(i),nz_pre(i),nt_next(i));
    %      B011(i) = B(nr_pre(i),nz_next(i),nt_next(i));
    %      B110(i) = B(nr_next(i),nz_next(i),nt_pre(i));
    %      B111(i) = B(nr_next(i),nz_next(i),nt_next(i));
    %     end
    % end
    if nr_next >= length(R) | nz_next >= length(Z)
        B_next = NaN;
    elseif nr_pre <= 0 | nz_pre <= 0
        B_next = NaN;
    else

        B000 = B(nr_pre,nz_pre,nt_pre);
        B100 = B(nr_next,nz_pre,nt_pre);
        B010 = B(nr_pre,nz_next,nt_pre);
        B001 = B(nr_pre,nz_pre,nt_next);
        B101 = B(nr_next,nz_pre,nt_next);
        B011 = B(nr_pre,nz_next,nt_next);
        B110 = B(nr_next,nz_next,nt_pre);
        B111 = B(nr_next,nz_next,nt_next);

        B_next = B000.*(1-rd).*(1-zd).*(1-td)+B100.*(1-zd).*(1-td)+B010.*(1-rd).*(1-td)+B001.*(1-rd).*(1-zd).*td+B101.*rd.*(1-zd).*td+B011.*(1-rd).*zd.*td+B110.*rd.*zd.*(1-td)+B111.*rd.*zd.*td;
    end

else

    B_next = nan;
end



