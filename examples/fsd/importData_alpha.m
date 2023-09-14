function importData_singleLatitude(P, pth)
    % Read FSD power law exponent (alpha). Output data in a separate file for
    % each year.
    %
    % Input array is dimensioned as [Year, Day].

    if nargin == 1
        pth = './data/raw';
    end

    dirIn = fullfile(pth, 'alpha');
    M = readmatrix(fullfile(dirIn, 'alphas_all.csv'));

    nD = 1;  % data space dimension
    idxDay1 = find(P.dayLim(1) == M(1, :), 1);
    idxDay2 = find(P.dayLim(2) == M(1, :), 1);
    dstr = sprintf('d%i-%i', P.dayLim);  % string for day range
    xystr = 'x1-1_y1-1';  % dummy string for lat/lon range
    xydstr = [xystr '_' dstr];

    nYr = P.yrLim(2) - P.yrLim(1) + 1;
    idxYr1 = find(P.yrLim(1) == M(:, 1), 1);

    for iYr = 1 : nYr
        i = idxYr1 + iYr - 1;
        x = M(i, idxDay1 : idxDay2);
        dirOut = fullfile(pth, ...
                          int2str(M(i, 1)), ...
                          'alpha', ...
                          xydstr);
        if ~isdir(dirOut)
            mkdir(dirOut)
        end
        save(fullfile(dirOut, 'dataGrid.mat'), 'nD')
        save(fullfile(dirOut, 'dataX.mat'), 'x')
    end
end
