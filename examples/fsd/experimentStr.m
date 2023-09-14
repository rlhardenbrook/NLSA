function s = experimentStr(P)
    % Make string identifier for data analysis experiment.
    s = strjoin_e({P.var ...
                   sprintf('%i-%i', P.yrLim(1), P.yrLim(2)) ...
                   sprintf('d%i-%i', P.dayLim(1), P.dayLim(2)) ...
                   sprintf('emb%i', P.embWindow) ...
                   P.kernel}, ...
                   '_');
    if isfield(P, 'den')
        if P.den
            s = [s '_den'];
        end
    end
end
