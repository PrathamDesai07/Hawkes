// Problem 8 simulation: relations between the regret notions and the parameter-error notions for the F-preconditioned online Newton estimator.
// Per run, at checkpoints n: realised log-loss regret r (vs truth), pseudo-regret R, quadratic sum Qsum = sum v'DJ v, per-gap sandwich extremes,
// information-weighted average theta_bar and (theta_bar-theta*)'J(theta_bar-theta*), last-iterate D_n = ‖theta-theta*‖_F^2 and e_n = ‖.‖_J^2,
// hindsight regret (batch MLE by Newton on the stored gap data).
// Model/algorithm as code/p6_scaling.cs (alpha*_k = (1-eps) beta_k / K, scaled prior, clip).
// Build: csc /o+ /out:p8_link.exe p8_link.cs   Run: p8_link.exe K=1 beta=1 mu=1 eps=0.1 reps=64 n=20000 start=near
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Threading.Tasks;

static class P8
{
    static int dbg = 0;
    static readonly double[] gx = { -0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498, 0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363 };
    static readonly double[] gw = { 0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620, 0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763 };

    static double[,] Chol(double[,] A, int n)
    {
        var L = new double[n, n];
        for (int i = 0; i < n; i++) for (int j = 0; j <= i; j++)
            { double s = A[i, j]; for (int k = 0; k < j; k++) s -= L[i, k] * L[j, k]; if (i == j) L[i, i] = Math.Sqrt(Math.Max(s, 1e-300)); else L[i, j] = s / L[j, j]; }
        return L;
    }
    static double[] Solve(double[,] A, int n, double[] b)
    {
        var L = Chol(A, n); var y = new double[n]; var x = new double[n];
        for (int i = 0; i < n; i++) { double s = b[i]; for (int k = 0; k < i; k++) s -= L[i, k] * y[k]; y[i] = s / L[i, i]; }
        for (int i = n - 1; i >= 0; i--) { double s = y[i]; for (int k = i + 1; k < n; k++) s -= L[k, i] * x[k]; x[i] = s / L[i, i]; }
        return x;
    }
    static double Quad(double[,] A, double[] v, int n) { double s = 0; for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) s += v[i] * A[i, j] * v[j]; return s; }
    static double U(Random r) { double u; do { u = r.NextDouble(); } while (u <= 0); return u; }

    class Out { public double[] r, R, Qs, Dn, en, avgErr, resid, hind, wilks, rmin, rmax, gapLo, gapHi, ratioBound, sandViol; }

    static Out Run(int K, double[] beta, double mu_s, double eps, int n, int seed, string start, double rho0, double warmEv, int[] cps)
    {
        int d = K + 1; var rng = new Random(seed); var al_s = new double[K]; for (int k = 0; k < K; k++) al_s[k] = (1 - eps) * beta[k] / K;
        var ths = new double[d]; ths[0] = mu_s; for (int k = 0; k < K; k++) ths[k + 1] = al_s[k];
        double lam_bar = mu_s / eps; var F = new double[d, d]; F[0, 0] = rho0 / lam_bar; for (int k = 0; k < K; k++) F[k + 1, k + 1] = rho0 * lam_bar / (beta[k] * beta[k]);
        var th = (double[])ths.Clone();
        if (start == "far") { th[0] = 2 * mu_s; for (int k = 0; k < K; k++) th[k + 1] = 0.5 * beta[k] / K; }
        else if (start == "near") { th[0] = 1.3 * mu_s; for (int k = 0; k < K; k++) th[k + 1] = 0.8 * al_s[k]; }
        var S = new double[K];
        Func<double> gap = () =>
        {
            double lam0 = mu_s; for (int k = 0; k < K; k++) lam0 += al_s[k] * S[k]; double w = 0;
            while (true) { w += -Math.Log(U(rng)) / lam0; double l = mu_s; for (int k = 0; k < K; k++) l += al_s[k] * S[k] * Math.Exp(-beta[k] * w); if (rng.NextDouble() * lam0 <= l) return w; }
        };
        for (long i = 0; i < (long)warmEv; i++) { double w = gap(); for (int k = 0; k < K; k++) S[k] = S[k] * Math.Exp(-beta[k] * w) + 1; }
        var A = new double[n][]; var B = new double[n][];     // stored gap data for the hindsight MLE
        var Jm = new double[d, d]; var sumJth = new double[d];
        double rr = 0, RR = 0, Qs = 0, bmax = beta.Max(); int G = cps.Length, gi = 0;
        double rmin = double.MaxValue, rmax = 0, ratioMin = double.MaxValue, ratioMax = 0;   // visited lam_theta/lam* range
        var gapR = new List<double[]>();   // (R_i, Q_i, ratio_min_i, ratio_max_i)
        var o = new Out { r = new double[G], R = new double[G], Qs = new double[G], Dn = new double[G], en = new double[G], avgErr = new double[G], resid = new double[G], hind = new double[G], wilks = new double[G], rmin = new double[G], rmax = new double[G], gapLo = new double[G], gapHi = new double[G], ratioBound = new double[G], sandViol = new double[G] };
        var phi = new double[d]; var a = new double[d]; var b = new double[d]; var dF = new double[d, d]; var dJ = new double[d, d];
        double sumRQlo = double.MaxValue, sumRQhi = 0, ellTheta = 0; double sandV = 0, sandC = 0; double cpl = 1, cml = 1;
        var th0log = 0.0; double ellSum = 0;
        double sumLoss = 0;   // sum_i l_i(theta_i) for hindsight regret
        for (int i = 1; i <= n; i++)
        {
            double w = gap(); a[0] = 1; b[0] = w;
            for (int k = 0; k < K; k++) { double e = Math.Exp(-beta[k] * w); a[k + 1] = S[k] * e; b[k + 1] = S[k] * (1 - e) / beta[k]; }
            A[i - 1] = (double[])a.Clone(); B[i - 1] = (double[])b.Clone();
            Array.Clear(dF, 0, dF.Length); Array.Clear(dJ, 0, dJ.Length);
            int nsub = Math.Max(2, Math.Min(12, (int)Math.Ceiling(w * bmax / 1.5))); double h = w / nsub, lg = 0, Qi = 0, rminG = double.MaxValue, rmaxG = 0;
            var v = new double[d]; for (int p = 0; p < d; p++) v[p] = th[p] - ths[p];
            for (int j = 0; j < nsub; j++) for (int q = 0; q < 8; q++)
                {
                    double u = h * (j + 0.5 * (gx[q] + 1)), wt = 0.5 * h * gw[q]; phi[0] = 1; for (int k = 0; k < K; k++) phi[k + 1] = S[k] * Math.Exp(-beta[k] * u);
                    double lt = 0, ls = 0, pv = 0; for (int p = 0; p < d; p++) { lt += phi[p] * th[p]; ls += phi[p] * ths[p]; pv += phi[p] * v[p]; }
                    lg += wt * ls * Math.Log(lt / ls); Qi += wt * pv * pv / ls; double ratio = lt / ls; rminG = Math.Min(rminG, ratio); rmaxG = Math.Max(rmaxG, ratio);
                    for (int p = 0; p < d; p++) for (int r2 = 0; r2 <= p; r2++) { double vv = phi[p] * phi[r2] * wt; dF[p, r2] += vv / lt; dJ[p, r2] += vv / ls; }
                }
            double bt = 0, bs = 0, at = 0, as_ = 0; for (int p = 0; p < d; p++) { bt += b[p] * th[p]; bs += b[p] * ths[p]; at += a[p] * th[p]; as_ += a[p] * ths[p]; }
            double Ri = bt - bs - lg; double ri = (bt - bs) - Math.Log(at / as_);   // realised regret increment: l_i(theta_i) - l_i(theta*)
            { double ra = at / as_; rminG = Math.Min(rminG, ra); rmaxG = Math.Max(rmaxG, ra);
              double l0t = th[0], l0s = ths[0]; for (int k = 0; k < K; k++) { l0t += th[k + 1] * S[k]; l0s += ths[k + 1] * S[k]; } double r0 = l0t / l0s; rminG = Math.Min(rminG, r0); rmaxG = Math.Max(rmaxG, r0); }
            rr += ri; RR += Ri; Qs += Qi; sumLoss += bt - Math.Log(at);
            ratioMin = Math.Min(ratioMin, rminG); ratioMax = Math.Max(ratioMax, rmaxG);
            if (Qi > 1e-9) { gapR.Add(new double[] { Ri, Qi, rminG, rmaxG, w }); }
            for (int p = 0; p < d; p++) for (int r2 = 0; r2 <= p; r2++) { Jm[p, r2] += dJ[p, r2]; if (r2 != p) Jm[r2, p] = Jm[p, r2]; }
            for (int p = 0; p < d; p++) { double s = 0; for (int r2 = 0; r2 < d; r2++) s += (r2 <= p ? dJ[p, r2] : dJ[r2, p]) * th[r2]; sumJth[p] += s; }
            for (int p = 0; p < d; p++) for (int r2 = 0; r2 <= p; r2++) { F[p, r2] += dF[p, r2]; if (r2 != p) F[r2, p] = F[p, r2]; }
            var g = new double[d]; for (int p = 0; p < d; p++) g[p] = b[p] - a[p] / at;
            var step = Solve(F, d, g); for (int p = 0; p < d; p++) th[p] -= step[p];
            th[0] = Math.Max(th[0], 0.1 * mu_s); double sa = 0; for (int k = 0; k < K; k++) { th[k + 1] = Math.Max(th[k + 1], 0.05 * al_s[k]); sa += th[k + 1] / beta[k]; }
            if (sa > 1 - eps / 2) for (int k = 0; k < K; k++) th[k + 1] *= (1 - eps / 2) / sa;
            for (int k = 0; k < K; k++) S[k] = S[k] * Math.Exp(-beta[k] * w) + 1;
            if (gi < G && i == cps[gi])
            {
                o.r[gi] = rr; o.R[gi] = RR; o.Qs[gi] = Qs;
                var vv2 = new double[d]; for (int p = 0; p < d; p++) vv2[p] = th[p] - ths[p];
                o.Dn[gi] = Quad(F, vv2, d); o.en[gi] = Quad(Jm, vv2, d);
                var tb = Solve(Jm, d, sumJth); var vb = new double[d]; for (int p = 0; p < d; p++) vb[p] = tb[p] - ths[p];
                o.avgErr[gi] = Quad(Jm, vb, d);
                // identity (3.1): Qs_total(with theta_i used at each gap) should equal avgErr + sum (theta_i - tb)'dJ(theta_i - tb); we check Qs >= avgErr and report the residual of the decomposition via sumJth
                o.resid[gi] = Qs - o.avgErr[gi];
                // sandwich extremes using range of ratios visited so far
                double cp = Math.Max(1, ratioMax), cm = Math.Min(1, ratioMin); o.rmin[gi] = ratioMin; o.rmax[gi] = ratioMax;
                double lo = double.MaxValue, hi = 0; int viol = 0;
                foreach (var gg in gapR) { double q = gg[0] / (gg[1] / 2); lo = Math.Min(lo, q); hi = Math.Max(hi, q);
                    double cpi = Math.Max(1, gg[3]), cmi = Math.Min(1, gg[2]);
                    if (gg[0] < gg[1] / (2 * cpi) * (1 - 1e-6) || gg[0] > gg[1] / (2 * cmi * cmi) * (1 + 1e-6)) { viol++; if (dbg < 6) { dbg++; Console.Error.WriteLine("VIOL R=" + gg[0].ToString("E4") + " Q/2=" + (gg[1] / 2).ToString("E4") + " lo=" + (gg[1] / (2 * cpi * cpi)).ToString("E4") + " hi=" + (gg[1] / (2 * cmi * cmi)).ToString("E4") + " rmin=" + gg[2].ToString("F4") + " rmax=" + gg[3].ToString("F4") + " w=" + gg[4].ToString("F3")); } } }
                o.gapLo[gi] = lo; o.gapHi[gi] = hi; o.sandViol[gi] = viol; o.ratioBound[gi] = 2 * cp * RR;
                // hindsight MLE by Newton with backtracking on L(theta) = sum (b'theta - ln a'theta)
                var tm = (double[])th.Clone(); double Lcur = 0;
                Func<double[], double> Lfun = (t) => { double L2 = 0; for (int j = 0; j < i; j++) { double ap = 0, bp = 0; for (int p = 0; p < d; p++) { ap += A[j][p] * t[p]; bp += B[j][p] * t[p]; } if (ap <= 1e-12) return double.PositiveInfinity; L2 += bp - Math.Log(ap); } return L2; };
                Lcur = Lfun(tm);
                for (int it = 0; it < 40; it++)
                {
                    var gr = new double[d]; var H = new double[d, d];
                    for (int j = 0; j < i; j++) { double ap = 0; for (int p = 0; p < d; p++) ap += A[j][p] * tm[p]; for (int p = 0; p < d; p++) { gr[p] += B[j][p] - A[j][p] / ap; for (int r2 = 0; r2 < d; r2++) H[p, r2] += A[j][p] * A[j][r2] / (ap * ap); } }
                    for (int p = 0; p < d; p++) H[p, p] += 1e-12;
                    var st = Solve(H, d, gr); double tstep = 1; bool ok = false;
                    for (int bt2 = 0; bt2 < 30; bt2++) { var tn = new double[d]; for (int p = 0; p < d; p++) tn[p] = tm[p] - tstep * st[p]; double Ln = Lfun(tn); if (Ln <= Lcur) { tm = tn; Lcur = Ln; ok = true; break; } tstep /= 2; }
                    if (!ok) break; double gn = 0; for (int p = 0; p < d; p++) gn += gr[p] * st[p]; if (gn < 1e-12) break;
                }
                o.hind[gi] = sumLoss - Lcur;               // hindsight regret = sum l_i(theta_i) - L(theta_hat)
                o.wilks[gi] = Lfun(ths) - Lcur;            // L(theta*) - L(theta_hat) >= 0, should average d/2
                gi++;
            }
        }
        return o;
    }

    static void Main(string[] args)
    {
        var kv = new Dictionary<string, string>(); foreach (var s in args) { var p = s.Split('='); kv[p[0]] = p[1]; }
        Func<string, string, string> get = (k, dflt) => kv.ContainsKey(k) ? kv[k] : dflt; var ci = CultureInfo.InvariantCulture;
        int K = int.Parse(get("K", "1")); var beta = get("beta", K == 1 ? "1" : "1,2").Split(',').Select(x => double.Parse(x, ci)).ToArray();
        double mu = double.Parse(get("mu", "1"), ci), eps = double.Parse(get("eps", "0.1"), ci); int reps = int.Parse(get("reps", "32")), n = int.Parse(get("n", "20000"));
        string start = get("start", "near"); double rho0 = double.Parse(get("rho0", "1"), ci); int seed0 = int.Parse(get("seed", "8000")); int d = K + 1;
        var cps = new List<int>(); for (double c = 100; c <= n * 1.0001; c *= Math.Pow(10, 0.5)) cps.Add((int)Math.Round(c)); if (cps.Last() != n) cps.Add(n); var cp = cps.ToArray();
        var all = new Out[reps]; double warm = 20.0 / (eps * eps);
        Parallel.For(0, reps, new ParallelOptions { MaxDegreeOfParallelism = 8 }, r => { all[r] = Run(K, beta, mu, eps, n, seed0 + r, start, rho0, warm, cp); });
        Console.WriteLine("# K=" + K + " beta=" + string.Join(",", beta) + " mu=" + mu + " eps=" + eps + " start=" + start + " reps=" + reps + " (d=" + d + ")");
        Console.WriteLine("# n  E_r(realised)  se  E_R(pseudo)  se  E_r-E_R  se_diff | E_Qsum/2  R/(Q/2) | E_Dn  E_en  E_avgErr  E_avgErr/(2*cp*R)max | hind_regret  E_wilks | gapR/(Q/2)min max  sandwich_viol  ratio_range")
        ;
        for (int g = 0; g < cp.Length; g++)
        {
            Func<Func<Out, double>, double> mean = f => all.Average(f);
            Func<Func<Out, double>, double> se = f => { double m = all.Average(f); return Math.Sqrt(all.Sum(x => (f(x) - m) * (f(x) - m)) / Math.Max(1, reps - 1) / reps); };
            double mr = mean(x => x.r[g]), mR = mean(x => x.R[g]);
            double sd = se(x => x.r[g] - x.R[g]);
            double maxFrac = all.Max(x => x.avgErr[g] / x.ratioBound[g]);
            Console.WriteLine(String.Format(ci, "{0,6} {1,8:F3} {2,6:F3} {3,8:F3} {4,6:F3} {5,8:F3} {6,6:F3} | {7,8:F3} {8,6:F3} | {9,7:F3} {10,7:F3} {11,8:F3} {12,6:F3} | {13,8:F3} {14,6:F3} | {15,6:F3} {16,6:F3} {17,5:F0} [{18:F3},{19:F3}] minres={20:E2}",
                cp[g], mr, se(x => x.r[g]), mR, se(x => x.R[g]), mr - mR, sd, mean(x => x.Qs[g]) / 2, mR / (mean(x => x.Qs[g]) / 2),
                mean(x => x.Dn[g]), mean(x => x.en[g]), mean(x => x.avgErr[g]), maxFrac, mean(x => x.hind[g]), mean(x => x.wilks[g]),
                all.Min(x => x.gapLo[g]), all.Max(x => x.gapHi[g]), all.Sum(x => x.sandViol[g]), all.Min(x => x.rmin[g]), all.Max(x => x.rmax[g]), all.Min(x => x.resid[g])));
        }
    }
}
