clear
clc
warning('off')

function phi = binaryMixture(theta, setFrac, phase)
    x1 = theta(1);
    P = theta(2);
    x = [x1; 1 - x1];
    y = [setFrac; 1 - setFrac];
    T = 307.4;
    R = 8.3144621e-3;
    Tc = [305.3; 660];
    Pc = [4872; 2750];
    omega = [0.1; 0.313];
    Tr = T ./ Tc;
    ki = 0.37464 + (1.54226 .* omega) - (0.26992 .* omega.^2);
    alpha = (1 + (ki .* (1 - sqrt(Tr)))).^2;
    ai = (0.45724 .* ((R.^2 .* Tc.^2) ./ Pc)) .* alpha;
    bi = 0.07780 * ((R .* Tc) ./ Pc);
    aij = sqrt(ai * ai');
    if(phase == "liquid")
        cp = x;
    elseif(phase == "vapor")
        cp = y;
    end
    a = cp' * aij * cp;
    A = (a * P) / (R.^2 * T.^2);
    b = cp' * bi;
    B = (b * P) / (R * T);
    eqState = poly([- ((A * B) - B.^2 - B.^3), A - (3 * B.^2) - (2 * B), - (1 - B), 1], "Z", "coeff");
    eqRoots = roots(eqState);
    boolRealRoots = eqRoots == real(eqRoots);
    realRoots = [];
    for i = 1 : length(boolRealRoots)
        if(boolRealRoots(i) == %T)
            realRoots($ + 1) = real(eqRoots(i));
        end
    end
    if(phase == "liquid")
        Z = min(realRoots);
    elseif(phase == "vapor")
        Z = max(realRoots);
    end
    phi = exp(((bi ./ b) .* (Z - 1)) - log(Z - B) - ((A ./ (2 .* sqrt(2) .* B)) .* (((2 ./ a) * (aij * cp)) - (bi ./ b)) .* log((Z + ((1 + sqrt(2)) .* B)) ./ (Z + ((1 - sqrt(2)) .* B)))));
endfunction

function F = sys(theta, setFrac, prob)
    if(prob == "bubble")
        x1 = setFrac;
        y1 = theta(1);
        P = theta(2);
        phiL = binaryMixture(theta, setFrac, "liquid");
        phiV = binaryMixture(theta, setFrac, "vapor");
    elseif(prob == "dew")
        y1 = setFrac;
        x1 = theta(1);
        P = theta(2);
        phiL = binaryMixture(theta, setFrac, "liquid");
        phiV = binaryMixture(theta, setFrac, "vapor");
    end
    F = [(x1 .* phiL(1)) - (y1 .* phiV(1)); ((1 - x1) .* phiL(2)) - ((1 - y1) .* phiV(2))];
endfunction

function y = diffFirstOrder(theta, m, setFrac, prob)
    dim = length(theta);
    a = zeros(dim, 1);
    b = zeros(dim, 1);
    c = zeros(dim, 1);
    d = zeros(dim, 1);
    h = 1e-6;
    for i = 1 : dim
        if i == m
            a(i) = theta(i) + h;
            b(i) = theta(i) + (2 * h);
            c(i) = theta(i) - h;
            d(i) = theta(i) - (2 * h);
        else
            a(i) = theta(i);
            b(i) = theta(i);
            c(i) = theta(i);
            d(i) = theta(i);
        end
    end
    y = ((8 .* sys(a, setFrac, prob)) - sys(b, setFrac, prob) - (8 .* sys(c, setFrac, prob)) + sys(d, setFrac, prob)) ./ (12 * h);
endfunction

function J = jacobian(theta, setFrac, prob)
    dim = length(theta);
    J = [];
    for i = 1 : dim
        J = [J, diffFirstOrder(theta, i, setFrac, prob)];
    end
endfunction

function thetaf = newtonMultidimRoot(theta, setFrac, prob)
    dim = length(theta);
    err = 1;
    threshold = 1e-6;
    thetaf = theta;
    thetai = thetaf;
    while err > threshold
        thetaf = thetai - (0.1 .* (inv(jacobian(thetaf, setFrac, prob)) * sys(thetaf, setFrac, prob)));
        err = norm(thetaf - thetai);
        thetai = thetaf;
    end
endfunction

function plotVLE(guess, problem)
    cont = 1;
    n = 500;
    epsilon = 1e-6;
    P = guess(2);
    res = [];
    if problem == "dew"
        x1 = guess(1);
        y1 = linspace(0.9989, 0.99898, n)';
        for i = 1 : n
            tempRes = newtonMultidimRoot([x1; P], y1(i), problem);
            if i > 1
                nanTest = sum(isnan(tempRes) .* [1; 1]);
                infTest = sum(isinf(tempRes) .* [1; 1]);
                realTest = isreal(tempRes);
                if (norm(res(:, i - 1) - tempRes) < 10) & (nanTest == 0) & (infTest == 0) & (realTest == %T)
                    res = [res, tempRes];
                    x1 = res(1, i);
                    P = res(2, i);
                    cont = cont + 1;
                else
                    break;
                end
            else
                res = [res, tempRes];
                x1 = res(1, i);
                P = res(2, i);
                cont = cont + 1;
            end
            disp(cont);
        end
        plot(y1(1 : cont - 1), res(2, :)');
    elseif problem == "bubble"
        x1 = linspace(0.9989, 0.99898, n)';
        y1 = guess(1);
        for i = 1 : n
            tempRes = newtonMultidimRoot([y1; P], x1(i), problem);
            if i > 1
                nanTest = sum(isnan(tempRes) .* [1; 1]);
                infTest = sum(isinf(tempRes) .* [1; 1]);
                realTest = isreal(tempRes);
                if (norm(res(:, i - 1) - tempRes) < 10) & (nanTest == 0) & (infTest == 0) & (realTest == %T)
                    res = [res, tempRes];
                    y1 = res(1, i);
                    P = res(2, i);
                    cont = cont + 1;
                else
                    break;
                end
            else
                res = [res, tempRes];
                y1 = res(1, i);
                P = res(2, i);
                cont = cont + 1;
            end
            disp(cont);
        end
        plot(x1(1 : cont - 1), res(2, :)');
    end
endfunction

guess = [0.999; 5015];
problem = "dew";
plotVLE(guess, problem);
