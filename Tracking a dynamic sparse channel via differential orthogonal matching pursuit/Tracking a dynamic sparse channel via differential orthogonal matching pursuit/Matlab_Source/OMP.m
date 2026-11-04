function h_f = OMP(L, Phi_f, y_f, k_max)
    h_f = zeros(L, 1);
    k_f = 0;
    r = y_f;
    while (k_f<k_max)
        k_f = k_f + 1;
        error = Phi_f'*r;
        [val, pos] = sort(abs(error), 'descend');
        if (k_f == 1)
            d_set = pos(1);
        else
            d_set = union(d_set, pos(1));
        end
        h_f(d_set) = pinv(Phi_f(:,d_set))*y_f;
        r = y_f - Phi_f*h_f;
    end
end