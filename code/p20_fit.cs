// Problem 20: exponential-dictionary Hawkes fits to the Deribit BTC option-trade stream (events.bin: float64 seconds since 2024-11-01 00:00 UTC, one entry per distinct millisecond).
// Window = [day*86400 + h0*3600, day*86400 + h1*3600]; the state is initialised from the events of the preceding 10 minutes (no likelihood contribution).
//   lam(t) = mu + sum_k alpha_k S_k(t), dS_k = -beta_k S_k dt + dN;  n_hat = sum_k alpha_k / beta_k.   beta_k fixed (per second).
// modes:  select  : pooled fit on train days, held-out log-lik per event on validation days, for a list of dictionaries
//         daily   : per-window batch MLE (several dictionaries), n_hat with delta-method s.e., time-rescaling KS and lag-1 correlation
//         online  : forgetting online Newton (eta = 1/(eta_mult * N_window_events)) through the days; end-of-window theta and n_hat
//         rolling : batch MLE on rolling windows of `wlen` hours over a day range (intra-day tracking)
// Build: csc /o+ /out:p20_fit.exe p20_fit.cs
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

static class P20
{
    static double[] ev;
    static CultureInfo ci = CultureInfo.InvariantCulture;
    static int LowerBound(double t) { int lo = 0, hi = ev.Length; while (lo < hi) { int m = (lo + hi) / 2; if (ev[m] < t) lo = m + 1; else hi = m; } return lo; }

    class Win { public double t0, t1; public int n; public double[] w; public double[][] a, b; public double[] bEnd; public double[][] S0; }
    // build statistics (a_i, b_i) for events of [t0,t1) with state initialised from the preceding 600 s
    static Win Build(double t0, double t1, double[] beta)
    {
        int K = beta.Length, d = K + 1; int i0 = LowerBound(t0), i1 = LowerBound(t1), iw = LowerBound(t0 - 600);
        var S = new double[K]; double tp = t0 - 600; for (int i = iw; i < i0; i++) { double w = ev[i] - tp; for (int k = 0; k < K; k++) S[k] = S[k] * Math.Exp(-beta[k] * w) + 1; tp = ev[i]; }
        { double w = t0 - tp; for (int k = 0; k < K; k++) S[k] *= Math.Exp(-beta[k] * w); }
        var W = new Win { t0 = t0, t1 = t1, n = i1 - i0, w = new double[i1 - i0], a = new double[i1 - i0][], b = new double[i1 - i0][], S0 = new double[i1 - i0][] };
        double last = t0;
        for (int i = i0; i < i1; i++)
        {
            int j = i - i0; double w = ev[i] - last; W.w[j] = w; W.S0[j] = (double[])S.Clone(); var a = new double[d]; var bb = new double[d]; a[0] = 1; bb[0] = w;
            for (int k = 0; k < K; k++) { double e = Math.Exp(-beta[k] * w); a[k + 1] = S[k] * e; bb[k + 1] = S[k] * (1 - e) / beta[k]; S[k] = S[k] * e + 1; }
            W.a[j] = a; W.b[j] = bb; last = ev[i];
        }
        { double w = t1 - last; var be = new double[d]; be[0] = w; for (int k = 0; k < K; k++) be[k + 1] = S[k] * (1 - Math.Exp(-beta[k] * w)) / beta[k]; W.bEnd = be; }
        return W;
    }
    static double[] Solve(double[,] A, int n, double[] b)
    {
        var L = new double[n, n];
        for (int i = 0; i < n; i++) for (int j = 0; j <= i; j++) { double s = A[i, j]; for (int k = 0; k < j; k++) s -= L[i, k] * L[j, k]; if (i == j) L[i, i] = Math.Sqrt(Math.Max(s, 1e-300)); else L[i, j] = s / L[j, j]; }
        var y = new double[n]; var x = new double[n];
        for (int i = 0; i < n; i++) { double s = b[i]; for (int k = 0; k < i; k++) s -= L[i, k] * y[k]; y[i] = s / L[i, i]; }
        for (int i = n - 1; i >= 0; i--) { double s = y[i]; for (int k = i + 1; k < n; k++) s -= L[k, i] * x[k]; x[i] = s / L[i, i]; }
        return x;
    }
    static double[,] Inverse(double[,] A, int n) { var inv = new double[n, n]; for (int c = 0; c < n; c++) { var e = new double[n]; e[c] = 1; var x = Solve(A, n, e); for (int r = 0; r < n; r++) inv[r, c] = x[r]; } return inv; }
    static void Clip(double[] th, double[] beta) { int K = beta.Length; th[0] = Math.Max(th[0], 1e-6); double sa = 0; for (int k = 0; k < K; k++) { th[k + 1] = Math.Max(th[k + 1], 1e-8); sa += th[k + 1] / beta[k]; } if (sa > 0.999) for (int k = 0; k < K; k++) th[k + 1] *= 0.999 / sa; }
    // negative log-likelihood sum (b'theta - ln a'theta) + b_end'theta over the given windows
    static double Nll(List<Win> ws, double[] th)
    {
        int d = th.Length; double L = 0;
        foreach (var W in ws) { for (int j = 0; j < W.n; j++) { double ap = 0, bp = 0; for (int q = 0; q < d; q++) { ap += W.a[j][q] * th[q]; bp += W.b[j][q] * th[q]; } if (ap <= 1e-300) return double.PositiveInfinity; L += bp - Math.Log(ap); } for (int q = 0; q < d; q++) L += W.bEnd[q] * th[q]; }
        return L;
    }
    static int NEv(List<Win> ws) { int s = 0; foreach (var W in ws) s += W.n; return s; }
    // constrained Newton with backtracking; returns theta and observed information (Hessian of nll)
    static double[] Fit(List<Win> ws, double[] beta, double[] start, out double[,] H, out double nll, out int iters)
    {
        int K = beta.Length, d = K + 1; var tm = (double[])start.Clone(); Clip(tm, beta); nll = Nll(ws, tm); iters = 0; H = new double[d, d];
        for (int it = 0; it < 80; it++)
        {
            var gr = new double[d]; H = new double[d, d];
            foreach (var W in ws) { for (int j = 0; j < W.n; j++) { double ap = 0; for (int q = 0; q < d; q++) ap += W.a[j][q] * tm[q]; double inv = 1 / ap; for (int q = 0; q < d; q++) { gr[q] += W.b[j][q] - W.a[j][q] * inv; for (int r = 0; r <= q; r++) H[q, r] += W.a[j][q] * W.a[j][r] * inv * inv; } } for (int q = 0; q < d; q++) gr[q] += W.bEnd[q]; }
            for (int q = 0; q < d; q++) { for (int r = 0; r < q; r++) H[r, q] = H[q, r]; H[q, q] += 1e-9 * (1 + H[q, q]); }
            var st = Solve(H, d, gr); double ts = 1; bool ok = false; iters++;
            for (int bt = 0; bt < 40; bt++) { var tn = new double[d]; for (int q = 0; q < d; q++) tn[q] = tm[q] - ts * st[q]; Clip(tn, beta); double Ln = Nll(ws, tn); if (Ln < nll - 1e-10 * Math.Abs(nll)) { tm = tn; nll = Ln; ok = true; break; } ts /= 2; }
            if (!ok) break;
        }
        return tm;
    }
    // time-rescaling residuals: compensator increments over gaps; returns KS vs Exp(1) and lag-1 correlation
    static void Resid(Win W, double[] th, double[] beta, out double ks, out double rho1, out double meanE)
    {
        int K = beta.Length; var E = new double[W.n]; int d = K + 1;
        for (int j = 0; j < W.n; j++) { double s = 0; for (int q = 0; q < d; q++) s += W.b[j][q] * th[q]; E[j] = s; }
        var srt = (double[])E.Clone(); Array.Sort(srt); double D = 0; int n = srt.Length; for (int i = 0; i < n; i++) { double F = 1 - Math.Exp(-srt[i]); D = Math.Max(D, Math.Max(Math.Abs(F - (double)i / n), Math.Abs((double)(i + 1) / n - F))); }
        ks = D; double m = E.Average(); meanE = m; double num = 0, den = 0; for (int i = 0; i < n; i++) { den += (E[i] - m) * (E[i] - m); if (i > 0) num += (E[i] - m) * (E[i - 1] - m); } rho1 = num / den;
    }
    static double[] Beta(string s) { return s.Split(',').Select(x => double.Parse(x, ci)).ToArray(); }
    static double[] Start(List<Win> ws, double[] beta) { int n = NEv(ws); double T = ws.Sum(w => w.t1 - w.t0); double lb = n / T; var th = new double[beta.Length + 1]; th[0] = 0.5 * lb; for (int k = 0; k < beta.Length; k++) th[k + 1] = 0.3 * beta[k] / beta.Length; return th; }
    static double Nhat(double[] th, double[] beta) { double n = 0; for (int k = 0; k < beta.Length; k++) n += th[k + 1] / beta[k]; return n; }
    static double NhatSe(double[,] H, double[] beta) { int K = beta.Length, d = K + 1; var Hi = Inverse(H, d); var c = new double[d]; for (int k = 0; k < K; k++) c[k + 1] = 1 / beta[k]; double v = 0; for (int p = 0; p < d; p++) for (int q = 0; q < d; q++) v += c[p] * Hi[p, q] * c[q]; return Math.Sqrt(Math.Max(v, 0)); }

    static void Main(string[] args)
    {
        var kv = new Dictionary<string, string>(); foreach (var s in args) { int p = s.IndexOf('='); kv[s.Substring(0, p)] = s.Substring(p + 1); }
        Func<string, string, string> get = (k, dflt) => kv.ContainsKey(k) ? kv[k] : dflt;
        string path = get("data", @"p20_data\events.bin"); var bytes = File.ReadAllBytes(path); ev = new double[bytes.Length / 8]; Buffer.BlockCopy(bytes, 0, ev, 0, bytes.Length);
        string mode = get("mode", "daily"); double h0 = double.Parse(get("h0", "13"), ci), h1 = double.Parse(get("h1", "20"), ci);
        var dicts = get("dicts", "0.1,1,10").Split(';').Select(Beta).ToArray();
        Func<int, Win> day = null;
        if (mode == "select")
        {
            int tr0 = int.Parse(get("tr0", "0")), tr1 = int.Parse(get("tr1", "30")), va0 = int.Parse(get("va0", "30")), va1 = int.Parse(get("va1", "60"));
            Console.WriteLine("# select: train days [" + tr0 + "," + tr1 + "), validation days [" + va0 + "," + va1 + "), window " + h0 + "-" + h1 + " UTC");
            Console.WriteLine("# dictionary(beta per s)  d  | train nll/ev  val ll/ev (pooled train fit) | val ll/ev (daily refit, in-sample) | pooled theta | n_hat");
            foreach (var beta in dicts)
            {
                var tr = new List<Win>(); for (int dd = tr0; dd < tr1; dd++) tr.Add(Build(dd * 86400.0 + h0 * 3600, dd * 86400.0 + h1 * 3600, beta));
                var va = new List<Win>(); for (int dd = va0; dd < va1; dd++) va.Add(Build(dd * 86400.0 + h0 * 3600, dd * 86400.0 + h1 * 3600, beta));
                double[,] H; double nl; int it; var th = Fit(tr, beta, Start(tr, beta), out H, out nl, out it);
                double vl = -Nll(va, th) / NEv(va); double loc = 0; int nloc = 0;
                foreach (var W in va) { var one = new List<Win> { W }; double[,] H2; double nl2; int it2; var th2 = Fit(one, beta, th, out H2, out nl2, out it2); loc += -nl2; nloc += W.n; }
                Console.WriteLine(String.Format(ci, "{0,-14} {1} | {2,9:F4} {3,9:F4} | {4,9:F4} | {5} | {6:F3}", string.Join(",", beta), beta.Length + 1, nl / NEv(tr), vl, loc / nloc, string.Join(" ", th.Select(x => x.ToString("G4", ci))), Nhat(th, beta)));
            }
        }
        else if (mode == "daily")
        {
            int d0 = int.Parse(get("d0", "0")), d1 = int.Parse(get("d1", "365"));
            Console.WriteLine("# daily: window " + h0 + "-" + h1 + " UTC; columns: day dict nev | nhat se_nhat | mu alpha_k... | nll/ev KS rho1 meanE");
            for (int dd = d0; dd < d1; dd++)
                foreach (var beta in dicts)
                {
                    var W = Build(dd * 86400.0 + h0 * 3600, dd * 86400.0 + h1 * 3600, beta); if (W.n < 500) continue; var one = new List<Win> { W };
                    double[,] H; double nl; int it; var th = Fit(one, beta, Start(one, beta), out H, out nl, out it); double ks, rho, me; Resid(W, th, beta, out ks, out rho, out me);
                    Console.WriteLine(String.Format(ci, "{0} {1} {2} | {3:F4} {4:F4} | {5} | {6:F4} {7:F4} {8:F3} {9:F3}", dd, string.Join(",", beta), W.n, Nhat(th, beta), NhatSe(H, beta), string.Join(" ", th.Select(x => x.ToString("G5", ci))), nl / W.n, ks, rho, me));
                }
        }
        else if (mode == "online")
        {
            var beta = dicts[0]; int K = beta.Length, d = K + 1; int d0 = int.Parse(get("d0", "0")), d1 = int.Parse(get("d1", "365")); double etaMult = double.Parse(get("etamult", "5"), ci); double rho0 = double.Parse(get("rho0", "1"), ci);
            Console.WriteLine("# online forgetting Newton eta = 1/(etamult * events per window of the previous day), etamult=" + etaMult + "; columns: day nev | nhat | mu alpha_k...");
            double[] th = null; double[,] F = null; double nPrev = 0;
            for (int dd = d0; dd < d1; dd++)
            {
                var W = Build(dd * 86400.0 + h0 * 3600, dd * 86400.0 + h1 * 3600, beta); if (W.n < 500) continue;
                if (th == null) { th = Start(new List<Win> { W }, beta); double lb = W.n / (W.t1 - W.t0); F = new double[d, d]; F[0, 0] = rho0 / lb; for (int k = 0; k < K; k++) F[k + 1, k + 1] = rho0 * lb / (beta[k] * beta[k]); nPrev = W.n; }
                double eta = etaMult > 0 ? 1.0 / (etaMult * nPrev) : 0.0; var g = new double[d]; var dF = new double[d, d];
                for (int j = 0; j < W.n; j++)
                {
                    // Fisher increment by 8-point Gauss-Legendre on the gap
                    Array.Clear(dF, 0, dF.Length); double w = W.w[j]; int nsub = Math.Max(2, Math.Min(12, (int)Math.Ceiling(w * beta.Max() / 1.5))); double h = w / nsub; var phi = new double[d];
                    double[] gx = { -0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498, 0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363 };
                    double[] gw = { 0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620, 0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763 };
                    for (int s = 0; s < nsub; s++) for (int q = 0; q < 8; q++) { double u = h * (s + 0.5 * (gx[q] + 1)), wt = 0.5 * h * gw[q]; phi[0] = 1; for (int k = 0; k < K; k++) phi[k + 1] = W.S0[j][k] * Math.Exp(-beta[k] * u); double lt = 0; for (int p = 0; p < d; p++) lt += phi[p] * th[p]; for (int p = 0; p < d; p++) for (int r = 0; r <= p; r++) dF[p, r] += wt * phi[p] * phi[r] / lt; }
                    double at = 0; for (int p = 0; p < d; p++) at += W.a[j][p] * th[p];
                    for (int p = 0; p < d; p++) for (int r = 0; r <= p; r++) { F[p, r] = (1 - eta) * F[p, r] + dF[p, r]; if (r != p) F[r, p] = F[p, r]; }
                    for (int p = 0; p < d; p++) g[p] = W.b[j][p] - W.a[j][p] / at;
                    var st = Solve(F, d, g); for (int p = 0; p < d; p++) th[p] -= st[p]; Clip(th, beta);
                }
                nPrev = W.n; Console.WriteLine(String.Format(ci, "{0} {1} | {2:F4} | {3}", dd, W.n, Nhat(th, beta), string.Join(" ", th.Select(x => x.ToString("G5", ci)))));
            }
        }
        else if (mode == "rolling")
        {
            var beta = dicts[0]; int dA = int.Parse(get("d0", "340")), dB = int.Parse(get("d1", "346")); double wlen = double.Parse(get("wlen", "1"), ci), step = double.Parse(get("step", "1"), ci);
            Console.WriteLine("# rolling " + wlen + " h windows, step " + step + " h; columns: t_hours_from_day0 nev | nhat se | mu alpha_k...");
            for (double t = dA * 24.0; t + wlen <= dB * 24.0 + 1e-9; t += step)
            {
                var W = Build(t * 3600, (t + wlen) * 3600, beta); if (W.n < 200) continue; var one = new List<Win> { W }; double[,] H; double nl; int it; var th = Fit(one, beta, Start(one, beta), out H, out nl, out it);
                Console.WriteLine(String.Format(ci, "{0:F1} {1} | {2:F4} {3:F4} | {4}", t, W.n, Nhat(th, beta), NhatSe(H, beta), string.Join(" ", th.Select(x => x.ToString("G5", ci)))));
            }
        }
    }
}
