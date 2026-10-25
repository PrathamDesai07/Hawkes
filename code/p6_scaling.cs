// Problem 6 simulation: time-resolved expected regret of the F-preconditioned online Newton estimator for a K-exponential Hawkes process.
// Model: lam = mu + sum_k alpha_k S_k, dS_k = -beta_k S_k dt + dN, alpha*_k = (1-eps) beta_k / K (so branching ratio n = 1-eps), mean intensity mu/eps.
// Algorithm (same as code/p5_k2_regret.py with a clip instead of the F-projection): F += int phi phi'/lam_theta, theta <- clip(theta - F^-1 g).
// Output per checkpoint time T: E R, s.e., E ln det F, E ln det F0, mean oracle information J/T (upper triangle), mean events.
// Build: csc /o+ /out:p6_scaling.exe p6_scaling.cs     Run:  p6_scaling.exe K=1 beta=1 mu=1 eps=0.1 reps=24 tmax=2000 start=far prior=scaled rho0=1
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Threading.Tasks;

static class P6
{
    static readonly double[] gx = { -0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498, 0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363 };
    static readonly double[] gw = { 0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620, 0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763 };

    static double SolveLogDet(double[,] A, int n, out double[,] L)
    {
        L = new double[n, n]; double ld = 0;
        for (int i = 0; i < n; i++)
            for (int j = 0; j <= i; j++)
            {
                double s = A[i, j];
                for (int k = 0; k < j; k++) s -= L[i, k] * L[j, k];
                if (i == j) { if (s <= 0) s = 1e-300; L[i, i] = Math.Sqrt(s); ld += 2 * Math.Log(L[i, i]); }
                else L[i, j] = s / L[j, j];
            }
        return ld;
    }
    static double[] CholSolve(double[,] L, int n, double[] b)
    {
        var y = new double[n]; var x = new double[n];
        for (int i = 0; i < n; i++) { double s = b[i]; for (int k = 0; k < i; k++) s -= L[i, k] * y[k]; y[i] = s / L[i, i]; }
        for (int i = n - 1; i >= 0; i--) { double s = y[i]; for (int k = i + 1; k < n; k++) s -= L[k, i] * x[k]; x[i] = s / L[i, i]; }
        return x;
    }
    static double Uniform(Random r) { double u; do { u = r.NextDouble(); } while (u <= 0); return u; }

    class Res { public double[] R, LD, T, N, DEF, TRS, CNT, GS, XS; public double[][] J; public double Texit; }

    static Res Run(int K, double[] beta, double mu_s, double eps, double tmax, int seed, string start, string prior, double rho0, double warmEv, double[] grid, double delta)
    {
        int d = K + 1; var rng = new Random(seed);
        var al_s = new double[K]; for (int k = 0; k < K; k++) al_s[k] = (1 - eps) * beta[k] / K;
        var ths = new double[d]; ths[0] = mu_s; for (int k = 0; k < K; k++) ths[k + 1] = al_s[k];
        double lam_bar = mu_s / eps;
        var F = new double[d, d];
        if (prior == "scaled") { F[0, 0] = rho0 / lam_bar; for (int k = 0; k < K; k++) F[k + 1, k + 1] = rho0 * lam_bar / (beta[k] * beta[k]); }
        else { F[0, 0] = rho0; for (int k = 0; k < K; k++) F[k + 1, k + 1] = rho0 / (beta[k] * beta[k]); }
        double[,] L0; double ld0 = SolveLogDet(F, d, out L0);
        var th = (double[])ths.Clone();
        if (start == "far") { th[0] = 2 * mu_s; for (int k = 0; k < K; k++) th[k + 1] = 0.5 * beta[k] / K; }
        else if (start == "near") { th[0] = 1.3 * mu_s; for (int k = 0; k < K; k++) th[k + 1] = 0.8 * al_s[k]; }
        var S = new double[K]; double Tt = 0, R = 0; int n = 0;
        var Jm = new double[d, d];
        Func<double> gap = () =>
        {
            double lam0 = mu_s; for (int k = 0; k < K; k++) lam0 += al_s[k] * S[k];
            double w = 0;
            while (true)
            {
                w += -Math.Log(Uniform(rng)) / lam0;
                double l = mu_s; for (int k = 0; k < K; k++) l += al_s[k] * S[k] * Math.Exp(-beta[k] * w);
                if (rng.NextDouble() * lam0 <= l) return w;
            }
        };
        for (long i = 0; i < (long)warmEv; i++) { double w = gap(); for (int k = 0; k < K; k++) S[k] = S[k] * Math.Exp(-beta[k] * w) + 1; }
        int G = grid.Length, gi = 0;
        var res = new Res { R = new double[G], LD = new double[G], T = new double[G], N = new double[G], DEF = new double[G], TRS = new double[G], CNT = new double[G], GS = new double[G], XS = new double[G], J = new double[G][] };
        double defect = 0, trsum = 0, cnt = 0, texit = 0, gsum = 0, xsum = 0;
        double bmax = beta.Max();
        var phi = new double[d]; var a = new double[d]; var b = new double[d]; var g = new double[d];
        var dF = new double[d, d]; var dJ = new double[d, d];
        double lo_mu = 0.1 * mu_s; double cap = 1 - eps / 2;
        while (gi < G)
        {
            double w = gap(); n++;
            a[0] = 1; b[0] = w;
            for (int k = 0; k < K; k++) { double e = Math.Exp(-beta[k] * w); a[k + 1] = S[k] * e; b[k + 1] = S[k] * (1 - e) / beta[k]; }
            Array.Clear(dF, 0, dF.Length); Array.Clear(dJ, 0, dJ.Length);
            int nsub = Math.Max(2, Math.Min(12, (int)Math.Ceiling(w * bmax / 1.5))); double h = w / nsub, lg = 0;
            for (int j = 0; j < nsub; j++)
                for (int q = 0; q < 8; q++)
                {
                    double u = h * (j + 0.5 * (gx[q] + 1)), wt = 0.5 * h * gw[q];
                    phi[0] = 1; for (int k = 0; k < K; k++) phi[k + 1] = S[k] * Math.Exp(-beta[k] * u);
                    double lt = 0, ls = 0; for (int p = 0; p < d; p++) { lt += phi[p] * th[p]; ls += phi[p] * ths[p]; }
                    lg += wt * ls * Math.Log(lt / ls);
                    for (int p = 0; p < d; p++) for (int r2 = 0; r2 <= p; r2++) { double v = phi[p] * phi[r2] * wt; dF[p, r2] += v / lt; dJ[p, r2] += v / ls; }
                }
            double bt = 0, bs = 0, at = 0; for (int p = 0; p < d; p++) { bt += b[p] * th[p]; bs += b[p] * ths[p]; at += a[p] * th[p]; }
            R += bt - bs - lg;
            double ldOld; double[,] Lo; ldOld = SolveLogDet(F, d, out Lo); double trX = 0;
            for (int c2 = 0; c2 < d; c2++) { var col = new double[d]; for (int p = 0; p < d; p++) col[p] = p >= c2 ? dF[p, c2] : dF[c2, p]; var sol = CholSolve(Lo, d, col); trX += sol[c2]; }
            { var gg = new double[d]; for (int p = 0; p < d; p++) gg[p] = b[p] - a[p] / at; var sg = CholSolve(Lo, d, gg); double gq = 0; for (int p = 0; p < d; p++) gq += gg[p] * sg[p]; gsum += gq;
              double xq = 0; for (int p = 0; p < d; p++) for (int r2 = 0; r2 < d; r2++) { double dfv = p >= r2 ? dF[p, r2] : dF[r2, p]; xq += (th[p] - ths[p]) * dfv * (th[r2] - ths[r2]); } xsum += xq; }
            for (int p = 0; p < d; p++) for (int r2 = 0; r2 <= p; r2++) { F[p, r2] += dF[p, r2]; if (r2 != p) F[r2, p] = F[p, r2]; Jm[p, r2] += dJ[p, r2]; }
            for (int p = 0; p < d; p++) g[p] = b[p] - a[p] / at;
            double[,] L; SolveLogDet(F, d, out L); var step = CholSolve(L, d, g);
            for (int p = 0; p < d; p++) th[p] -= step[p];
            th[0] = Math.Max(th[0], lo_mu);
            double sa = 0; for (int k = 0; k < K; k++) { th[k + 1] = Math.Max(th[k + 1], 0.05 * al_s[k]); sa += th[k + 1] / beta[k]; }
            if (sa > cap) for (int k = 0; k < K; k++) th[k + 1] *= cap / sa;
            for (int k = 0; k < K; k++) S[k] = S[k] * Math.Exp(-beta[k] * w) + 1;
            Tt += w;
            { double[,] Ln; double ldNew = SolveLogDet(F, d, out Ln); defect += trX - (ldNew - ldOld); trsum += trX; if (trX > 1) cnt++; }
            { bool out0 = Math.Abs(th[0] / mu_s - 1) > delta; for (int k = 0; k < K; k++) if (Math.Abs(th[k + 1] / al_s[k] - 1) > delta) out0 = true; if (out0) texit = Tt; }
            while (gi < G && Tt >= grid[gi])
            {
                double[,] L2; double ld = SolveLogDet(F, d, out L2);
                res.R[gi] = R; res.LD[gi] = ld - ld0; res.T[gi] = Tt; res.N[gi] = n; res.DEF[gi] = defect; res.TRS[gi] = trsum; res.CNT[gi] = cnt; res.GS[gi] = gsum; res.XS[gi] = xsum;
                var jv = new double[d * (d + 1) / 2]; int c = 0; for (int p = 0; p < d; p++) for (int r2 = 0; r2 <= p; r2++) jv[c++] = Jm[p, r2];
                res.J[gi] = jv; gi++;
            }
        }
        res.Texit = texit; return res;
    }

    static double[] Eig2or3(double[,] A, int d)
    {   // Jacobi eigenvalues for small symmetric matrices
        var M = (double[,])A.Clone();
        for (int sweep = 0; sweep < 60; sweep++)
        {
            double off = 0; for (int i = 0; i < d; i++) for (int j = 0; j < i; j++) off += M[i, j] * M[i, j];
            if (off < 1e-24) break;
            for (int p = 0; p < d; p++) for (int q = p + 1; q < d; q++)
                {
                    if (Math.Abs(M[p, q]) < 1e-300) continue;
                    double th = 0.5 * Math.Atan2(2 * M[p, q], M[q, q] - M[p, p]), c = Math.Cos(th), s = Math.Sin(th);
                    for (int k = 0; k < d; k++) { double kp = M[k, p], kq = M[k, q]; M[k, p] = c * kp - s * kq; M[k, q] = s * kp + c * kq; }
                    for (int k = 0; k < d; k++) { double pk = M[p, k], qk = M[q, k]; M[p, k] = c * pk - s * qk; M[q, k] = s * pk + c * qk; }
                }
        }
        var ev = new double[d]; for (int i = 0; i < d; i++) ev[i] = M[i, i]; Array.Sort(ev); return ev;
    }

    static void Main(string[] args)
    {
        var kv = new Dictionary<string, string>(); foreach (var s in args) { var p = s.Split('='); kv[p[0]] = p[1]; }
        Func<string, string, string> get = (k, dflt) => kv.ContainsKey(k) ? kv[k] : dflt;
        var ci = CultureInfo.InvariantCulture;
        int K = int.Parse(get("K", "1")); var beta = get("beta", K == 1 ? "1" : "1,2").Split(',').Select(x => double.Parse(x, ci)).ToArray();
        double mu = double.Parse(get("mu", "1"), ci), eps = double.Parse(get("eps", "0.1"), ci), tmax = double.Parse(get("tmax", "1000"), ci), tmin = double.Parse(get("tmin", "1"), ci);
        int reps = int.Parse(get("reps", "16")); string start = get("start", "far"), prior = get("prior", "scaled"); double rho0 = double.Parse(get("rho0", "1"), ci);
        double warm = double.Parse(get("warm", "0"), ci) ; if (warm == 0) warm = 20.0 / (eps * eps);
        int seed0 = int.Parse(get("seed", "1000")); double delta = double.Parse(get("delta", "0.25"), ci); int d = K + 1;
        var grid = new List<double>(); for (double t = tmin; t <= tmax * 1.0001; t *= Math.Pow(10, 0.25)) grid.Add(t);
        var G = grid.ToArray();
        var all = new Res[reps];
        Parallel.For(0, reps, new ParallelOptions { MaxDegreeOfParallelism = 8 }, r => { all[r] = Run(K, beta, mu, eps, tmax, seed0 + r, start, prior, rho0, warm, G, delta); });
        double k0 = 2 * mu / beta[0];
        Console.WriteLine("# K=" + K + " beta=" + string.Join(",", beta) + " mu=" + mu + " eps=" + eps + " start=" + start + " prior=" + prior + " rho0=" + rho0 + " reps=" + reps + " warm_events=" + warm.ToString("F0") + " (d/2=" + (d / 2.0) + ")");
        Console.WriteLine("# T  E_R  se  R-(d/2)lnT  E_lndetF  lam_min(J/T)  lam_max(J/T)  nEvents  R-(d/2)lnN  tier_pred  E_defect  E_trace  E_count(x>1)  E_G  E_X");
        for (int gi = 0; gi < G.Length; gi++)
        {
            double mR = all.Average(x => x.R[gi]), mLD = all.Average(x => x.LD[gi]), mT = all.Average(x => x.T[gi]), mN = all.Average(x => x.N[gi]);
            double se = Math.Sqrt(all.Sum(x => (x.R[gi] - mR) * (x.R[gi] - mR)) / Math.Max(1, reps - 1) / reps);
            var Jb = new double[d, d]; int c = 0; for (int p = 0; p < d; p++) for (int q = 0; q <= p; q++) { double v = all.Average(x => x.J[gi][c]) / mT; Jb[p, q] = v; Jb[q, p] = v; c++; }
            var ev = Eig2or3(Jb, d);
            // tier-resolved prediction 0.5*sum ln(1 + T lam_j / f_j) with f_j = eigenvalue of the prior in the scaled-prior case along the matching axis (approximation: use f=rho0*scale of the diagonal prior entries sorted)
            double lam_bar = mu / eps; var f = new List<double>(); if (prior == "scaled") { f.Add(rho0 / lam_bar); for (int k = 0; k < K; k++) f.Add(rho0 * lam_bar / (beta[k] * beta[k])); } else { f.Add(rho0); for (int k = 0; k < K; k++) f.Add(rho0 / (beta[k] * beta[k])); }
            f.Sort(); double pred = 0; for (int j = 0; j < d; j++) pred += 0.5 * Math.Log(1 + mT * ev[j] / f[j]);
            Console.WriteLine(String.Format(ci, "{0,10:F1} {1,8:F3} {2,6:F3} {3,8:F3} {4,8:F3} {5,10:E3} {6,10:E3} {7,9:F0} {8,8:F3} {9,8:F3} {10,9:F3} {11,9:F2} {12,8:F2} {13,8:F2} {14,8:F2}", mT, mR, se, mR - d / 2.0 * Math.Log(mT), mLD * 0.5, ev[0], ev[d - 1], mN, mR - d / 2.0 * Math.Log(mN), pred, all.Average(x => x.DEF[gi]), all.Average(x => x.TRS[gi]), all.Average(x => x.CNT[gi]), all.Average(x => x.GS[gi]), all.Average(x => x.XS[gi])));
        }
        var te = all.Select(x => x.Texit).OrderBy(x => x).ToArray();
        Console.WriteLine("# Texit(delta=" + delta + ") mean=" + te.Average().ToString("F1") + " median=" + te[te.Length / 2].ToString("F1") + " p90=" + te[(int)(0.9 * (te.Length - 1))].ToString("F1"));
    }
}
