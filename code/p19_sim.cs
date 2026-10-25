// Problem 19: simulation study of the online Newton estimator with Ogata thinning and KNOWN (mu, alpha, beta).
// Emits RAW per-replicate rows at geometric time checkpoints so that slopes / confidence intervals are computed afterwards (code/p19_analyze.py).
//   row:  rep  T  n  R  D  e  v_0..v_{d-1}  J_0..J_{d(d+1)/2-1}      (R: pseudo-regret sum_i int lam* psi(lam_theta/lam*) du up to the last event before T;
//         D = (theta-theta*)' F (theta-theta*);  e = (theta-theta*)' J (theta-theta*), J = oracle information at theta*;  v = theta - theta*;  J upper triangle accumulated)
//   plus one line per replicate:  #TEXIT rep Texit_delta   (last time the iterate was outside R(delta): |mu/mu*-1|<=delta and |alpha_k/alpha*_k-1|<=delta)
// Model: lam = mu + sum_k alpha_k S_k, dS_k = -beta_k S_k dt + dN, alpha*_k = (1-eps) beta_k/K, stationary start (warm-up 20/eps^2 events), prior F_0 = rho0 diag(1/lam_bar, lam_bar/beta_k^2).
// Build: csc /o+ /out:p19_sim.exe p19_sim.cs   Run: p19_sim.exe K=1 beta=1 mu=1 eps=0.1 reps=200 tmax=30000 start=near > out.csv
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

static class P19
{
    static readonly double[] gx = { -0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498, 0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363 };
    static readonly double[] gw = { 0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620, 0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763 };
    static double U(Random r) { double u; do { u = r.NextDouble(); } while (u <= 0); return u; }
    static double[] SolveSPD(double[,] A, int n, double[] b)
    {
        var L = new double[n, n];
        for (int i = 0; i < n; i++) for (int j = 0; j <= i; j++) { double s = A[i, j]; for (int k = 0; k < j; k++) s -= L[i, k] * L[j, k]; if (i == j) L[i, i] = Math.Sqrt(Math.Max(s, 1e-300)); else L[i, j] = s / L[j, j]; }
        var y = new double[n]; var x = new double[n];
        for (int i = 0; i < n; i++) { double s = b[i]; for (int k = 0; k < i; k++) s -= L[i, k] * y[k]; y[i] = s / L[i, i]; }
        for (int i = n - 1; i >= 0; i--) { double s = y[i]; for (int k = i + 1; k < n; k++) s -= L[k, i] * x[k]; x[i] = s / L[i, i]; }
        return x;
    }
    static double Quad(double[,] A, double[] v, int n) { double s = 0; for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) s += v[i] * A[i, j] * v[j]; return s; }

    static string Run(int rep, int K, double[] beta, double mu_s, double eps, double tmax, int seed, string start, double rho0, double[] grid, out double texit, double delta)
    {
        int d = K + 1; var rng = new Random(seed); var al_s = new double[K]; for (int k = 0; k < K; k++) al_s[k] = (1 - eps) * beta[k] / K;
        var ths = new double[d]; ths[0] = mu_s; for (int k = 0; k < K; k++) ths[k + 1] = al_s[k];
        double lam_bar = mu_s / eps; var F = new double[d, d]; F[0, 0] = rho0 / lam_bar; for (int k = 0; k < K; k++) F[k + 1, k + 1] = rho0 * lam_bar / (beta[k] * beta[k]);
        var th = (double[])ths.Clone();
        if (start == "far") { th[0] = 2 * mu_s; for (int k = 0; k < K; k++) th[k + 1] = 0.5 * beta[k] / K; }
        else if (start == "near") { th[0] = 1.3 * mu_s; for (int k = 0; k < K; k++) th[k + 1] = 0.8 * al_s[k]; }
        var S = new double[K];
        Func<double> gap = () => { double lam0 = mu_s; for (int k = 0; k < K; k++) lam0 += al_s[k] * S[k]; double w = 0;
            while (true) { w += -Math.Log(U(rng)) / lam0; double l = mu_s; for (int k = 0; k < K; k++) l += al_s[k] * S[k] * Math.Exp(-beta[k] * w); if (rng.NextDouble() * lam0 <= l) return w; } };
        int warm = (int)(20.0 / (eps * eps)); for (int i = 0; i < warm; i++) { double w = gap(); for (int k = 0; k < K; k++) S[k] = S[k] * Math.Exp(-beta[k] * w) + 1; }
        var Jm = new double[d, d]; double R = 0, Tt = 0, bmax = beta.Max(); int n = 0, gi = 0, G = grid.Length; texit = 0;
        var phi = new double[d]; var a = new double[d]; var b = new double[d]; var dF = new double[d, d]; var dJ = new double[d, d]; var sb = new StringBuilder();
        double ca = 1 - eps / 2;
        while (gi < G)
        {
            double w = gap(); n++; a[0] = 1; b[0] = w;
            for (int k = 0; k < K; k++) { double e = Math.Exp(-beta[k] * w); a[k + 1] = S[k] * e; b[k + 1] = S[k] * (1 - e) / beta[k]; }
            Array.Clear(dF, 0, dF.Length); Array.Clear(dJ, 0, dJ.Length);
            int nsub = Math.Max(2, Math.Min(12, (int)Math.Ceiling(w * bmax / 1.5))); double h = w / nsub, lg = 0;
            for (int j = 0; j < nsub; j++) for (int q = 0; q < 8; q++)
                {
                    double u = h * (j + 0.5 * (gx[q] + 1)), wt = 0.5 * h * gw[q]; phi[0] = 1; for (int k = 0; k < K; k++) phi[k + 1] = S[k] * Math.Exp(-beta[k] * u);
                    double lt = 0, ls = 0; for (int p = 0; p < d; p++) { lt += phi[p] * th[p]; ls += phi[p] * ths[p]; }
                    lg += wt * ls * Math.Log(lt / ls);
                    for (int p = 0; p < d; p++) for (int r2 = 0; r2 <= p; r2++) { double vv = phi[p] * phi[r2] * wt; dF[p, r2] += vv / lt; dJ[p, r2] += vv / ls; }
                }
            double bt = 0, bs = 0, at = 0; for (int p = 0; p < d; p++) { bt += b[p] * th[p]; bs += b[p] * ths[p]; at += a[p] * th[p]; }
            R += bt - bs - lg;
            for (int p = 0; p < d; p++) for (int r2 = 0; r2 <= p; r2++) { F[p, r2] += dF[p, r2]; if (r2 != p) F[r2, p] = F[p, r2]; Jm[p, r2] += dJ[p, r2]; if (r2 != p) Jm[r2, p] = Jm[p, r2]; }
            var g = new double[d]; for (int p = 0; p < d; p++) g[p] = b[p] - a[p] / at;
            var step = SolveSPD(F, d, g); for (int p = 0; p < d; p++) th[p] -= step[p];
            th[0] = Math.Max(th[0], 0.1 * mu_s); double sa = 0; for (int k = 0; k < K; k++) { th[k + 1] = Math.Max(th[k + 1], 0.05 * al_s[k]); sa += th[k + 1] / beta[k]; }
            if (sa > ca) for (int k = 0; k < K; k++) th[k + 1] *= ca / sa;
            for (int k = 0; k < K; k++) S[k] = S[k] * Math.Exp(-beta[k] * w) + 1; Tt += w;
            bool outside = Math.Abs(th[0] / mu_s - 1) > delta; for (int k = 0; k < K; k++) if (Math.Abs(th[k + 1] / al_s[k] - 1) > delta) outside = true; if (outside) texit = Tt;
            while (gi < G && Tt >= grid[gi])
            {
                var v = new double[d]; for (int p = 0; p < d; p++) v[p] = th[p] - ths[p];
                sb.Append(rep).Append(' ').Append(Tt.ToString("R", CultureInfo.InvariantCulture)).Append(' ').Append(n).Append(' ').Append(R.ToString("R", CultureInfo.InvariantCulture)).Append(' ')
                  .Append(Quad(F, v, d).ToString("R", CultureInfo.InvariantCulture)).Append(' ').Append(Quad(Jm, v, d).ToString("R", CultureInfo.InvariantCulture));
                for (int p = 0; p < d; p++) sb.Append(' ').Append(v[p].ToString("R", CultureInfo.InvariantCulture));
                for (int p = 0; p < d; p++) for (int r2 = 0; r2 <= p; r2++) sb.Append(' ').Append(Jm[p, r2].ToString("R", CultureInfo.InvariantCulture));
                sb.Append('\n'); gi++;
            }
        }
        return sb.ToString();
    }

    static void Main(string[] args)
    {
        var kv = new Dictionary<string, string>(); foreach (var s in args) { var p = s.Split('='); kv[p[0]] = p[1]; }
        Func<string, string, string> get = (k, dflt) => kv.ContainsKey(k) ? kv[k] : dflt; var ci = CultureInfo.InvariantCulture;
        int K = int.Parse(get("K", "1")); var beta = get("beta", K == 1 ? "1" : "1,2").Split(',').Select(x => double.Parse(x, ci)).ToArray();
        double mu = double.Parse(get("mu", "1"), ci), eps = double.Parse(get("eps", "0.1"), ci), tmax = double.Parse(get("tmax", "10000"), ci), tmin = double.Parse(get("tmin", "1"), ci);
        int reps = int.Parse(get("reps", "100")), seed0 = int.Parse(get("seed", "19000")); string start = get("start", "near"); double rho0 = double.Parse(get("rho0", "1"), ci), delta = double.Parse(get("delta", "0.25"), ci);
        var grid = new List<double>(); for (double t = tmin; t <= tmax * 1.0001; t *= Math.Pow(10, 0.125)) grid.Add(t); var G = grid.ToArray();
        var outs = new string[reps]; var tex = new double[reps];
        Parallel.For(0, reps, new ParallelOptions { MaxDegreeOfParallelism = 8 }, r => { double te; outs[r] = Run(r, K, beta, mu, eps, tmax, seed0 + r, start, rho0, G, out te, delta); tex[r] = te; });
        var w = new System.IO.StreamWriter(Console.OpenStandardOutput());
        w.Write("# K=" + K + " beta=" + string.Join(",", beta) + " mu=" + mu.ToString(ci) + " eps=" + eps.ToString(ci) + " start=" + start + " rho0=" + rho0.ToString(ci) + " reps=" + reps + " delta=" + delta.ToString(ci) + " d=" + (K + 1) + "\n");
        w.Write("# rep T n R D e v_0.. J_upper_triangle(row-major lower)\n");
        for (int r = 0; r < reps; r++) w.Write(outs[r]);
        for (int r = 0; r < reps; r++) w.Write("#TEXIT " + r + " " + tex[r].ToString("R", ci) + "\n");
        w.Flush();
    }
}
