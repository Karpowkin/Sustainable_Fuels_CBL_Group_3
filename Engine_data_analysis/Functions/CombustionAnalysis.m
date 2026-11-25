function CC = CombustionAnalysis(Ca, p, V, SpS, Opt)
% aROHR analysis for one cycle. Returns CA10, CA50, CA90 and related quantities.

global Runiv
Ca = Ca(:); p = p(:); V = V(:);
dCa = Ca(2) - Ca(1);

% pegging
CC.pegOffset = 0;
mode = 'none';
if isfield(Opt,'PegMode'), mode = Opt.PegMode; end
switch mode
    case 'polytropic'
        i1 = abs(Ca - Opt.CaPoly(1)) <= 1;
        i2 = abs(Ca - Opt.CaPoly(2)) <= 1;
        p1 = mean(p(i1)); p2 = mean(p(i2));
        V1 = mean(V(i1)); V2 = mean(V(i2));
        k = Opt.kPoly;
        off = (p2*V2^k - p1*V1^k) / (V1^k - V2^k);
        p = p + off; CC.pegOffset = off;
    case 'ambient'
        iPeg = abs(Ca - Opt.CaPeg) <= 3;
        off = -mean(p(iPeg)) + Opt.pIntake;
        p = p + off; CC.pegOffset = off;
    case 'none'
        % no pegging
end

% filtering and derivatives
nWin  = 2*round(Opt.FiltWin/dCa/2) + 1;
pf    = smoothdata(p, 'sgolay', nWin);
dpdCa = gradient(pf, dCa);
dVdCa = gradient(V,  dCa);

% ratio of specific heats
iO2  = find(strcmp({SpS.Name}, 'O2'));
iN2  = find(strcmp({SpS.Name}, 'N2'));
Y    = [0.232 0.768];
Mair = 1/(Y(1)/SpS(iO2).Mass + Y(2)/SpS(iN2).Mass);
Rair = Runiv/Mair;
[~, iIVC] = min(abs(Ca - Opt.CaIVC));
m = pf(iIVC)*V(iIVC)/(Rair*Opt.TIVC);
T = pf.*V/(m*Rair);
if ischar(Opt.Gamma)
    Tc  = min(max(T, 200), 3000);
    cp  = Y(1)*CpNasa(Tc, SpS(iO2)) + Y(2)*CpNasa(Tc, SpS(iN2));
    gam = cp./(cp - Rair);
else
    gam = Opt.Gamma*ones(size(Ca));
end

% heat release
aROHR = gam./(gam-1).*pf.*dVdCa + 1./(gam-1).*V.*dpdCa;
iw    = Ca >= Opt.CaIVC & Ca <= Opt.CaEVO;
aHR   = nan(size(Ca));
aHR(iw) = cumtrapz(Ca(iw), aROHR(iw));

% mass fraction burned
CaW = Ca(iw);  aHRW = aHR(iw);  ROHRW = aROHR(iw);
iSearch = find(CaW >= Opt.CaSOI - 10);
[~, k] = max(ROHRW(iSearch));  iPeak = iSearch(k);
[~, iMin] = min(aHRW(iSearch(1):iPeak));  iMin = iSearch(1) - 1 + iMin;
[~, iMax] = max(aHRW(iMin:end));          iMax = iMin - 1 + iMax;
Qtot = aHRW(iMax) - aHRW(iMin);
xbW  = (aHRW - aHRW(iMin))/Qtot;
xb   = nan(size(Ca));  xb(iw) = xbW;
CaX  = @(X) crossing(CaW(iMin:iMax), xbW(iMin:iMax), X);

% results
CC.CA10          = CaX(0.10);
CC.CA50          = CaX(0.50);
CC.CA90          = CaX(0.90);
CC.BurnDuration  = CC.CA90 - CC.CA10;
CC.IgnitionDelay = CC.CA10 - Opt.CaSOI;
CC.Qtotal        = Qtot;
[CC.pMax, ip]    = max(pf);      CC.CApMax = Ca(ip);
CC.dpdCaMax      = max(dpdCa);
CC.aROHRMax      = ROHRW(iPeak); CC.CAaROHRMax = CaW(iPeak);
CC.Ca = Ca; CC.p = pf; CC.aROHR = aROHR; CC.aHR = aHR; CC.xb = xb; CC.gamma = gam; CC.T = T;
CC.Rair = Rair;
end

function CaX = crossing(Ca, xb, X)
k = find(xb >= X, 1);
if isempty(k) || k == 1
    CaX = NaN; return
end
CaX = Ca(k-1) + (X - xb(k-1))*(Ca(k) - Ca(k-1))/(xb(k) - xb(k-1));
end
