// Problem 21: baselines on identical simulated paths.  Methods (all use the same box clip and the same "near" start):
//   newton  : Fisher-preconditioned online Newton (ours, Problem 5), scaled prior rho0 = 1
//   batch   : full-sample MLE by Newton with backtracking, refit at the checkpoints (warm start = previous refit)      [no sequential regret]
//   hw_{c|s|i}<eta0> : Hall-Willett-type mirror descent with the log (multiplicative) mirror map on theta, step eta0 (c) / eta0/sqrt(i) (s) / eta0/i (i)
//   ls_r{1|2000|200} : Solo-type exponentially weighted least squares (forgetting rho = 1, 1-1/2000, 1-1/200 per event)
// mode=acc : raw rows  rep method n T R e v_0..v_{d-1}   (R = pseudo-regret of the sequential predictor, NaN for batch; e = (theta-theta*)'J_n(theta-theta*), J_n oracle)
//            plus a line  #I rep n T  J_upper_triangle   (oracle information for the normalisation)
// mode=time : sequential single-thread timing of one full pass per method (microseconds per event), batch refit time vs n.
// Build: csc /o+ /out:p21_baselines.exe p21_baselines.cs     Run: p21_baselines.exe mode=acc K=1 eps=0.1 reps=64 n=100000
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

static class P21
{
    static readonly double[] gx = { -0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498, 0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363 };
    static readonly double[] gw = { 0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620, 0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763 };
    static double U(Random r) { double u; do { u = r.NextDouble(); } while (u <= 0); return u; }
    static double[] Solve(double[,] A, int n, double[] b)
    {
        var L = new double[n, n];
        for (int i = 0; i < n; i++) for (int j = 0; j <= i; j++) { double s = A[i, j]; for (int k = 0; k < j; k++) s -= L[i, k] * L[j, k]; if (i == j) L[i, i] = Math.Sqrt(Math.Max(s, 1e-300)); else L[i, j] = s / L[j, j]; }
        var y = new double[n]; var x = new double[n];
        for (int i = 0; i < n; i++) { double s = b[i]; for (int k = 0; k < i; k++) s -= L[i, k] * y[k]; y[i] = s / L[i, i]; }
        for (int i = n - 1; i >= 0; i--) { double s = y[i]; for (int k = i + 1; k < n; k++) s -= L[k, i] * x[k]; x[i] = s / L[i, i]; }
        return x;
    }
    class Cfg { public int K, d; public double[] beta, ths, al_s, th0; public double mu_s, eps, lam_bar, rho0; public double[,] F0; public double lo_mu, cap; public double[] alMin; }
    class Path { public int n; public double[] w; public double[][] S; public double[][] a, b; public double T; }

    static Path Gen(Cfg c, int n, int seed)
    {
        var rng = new Random(seed); int K = c.K; var S = new double[K]; var p = new Path { n = n, w = new double[n], S = new double[n][], a = new double[n][], b = new double[n][] };
        Func<double> gap = () => { double lam0 = c.mu_s; for (int k = 0; k < K; k++) lam0 += c.al_s[k] * S[k]; double w = 0;
            while (true) { w += -Math.Log(U(rng)) / lam0; double l = c.mu_s; for (int k = 0; k < K; k++) l += c.al_s[k] * S[k] * Math.Exp(-c.beta[k] * w); if (rng.NextDouble() * lam0 <= l) return w; } };
        int warm = (int)(20.0 / (c.eps * c.eps)); for (int i = 0; i < warm; i++) { double w = gap(); for (int k = 0; k < K; k++) S[k] = S[k] * Math.Exp(-c.beta[k] * w) + 1; }
        double T = 0;
        for (int i = 0; i < n; i++)
        {
            double w = gap(); p.w[i] = w; p.S[i] = (double[])S.Clone(); var a = new double[c.d]; var b = new double[c.d]; a[0] = 1; b[0] = w;
            for (int k = 0; k < K; k++) { double e = Math.Exp(-c.beta[k] * w); a[k + 1] = S[k] * e; b[k + 1] = S[k] * (1 - e) / c.beta[k]; S[k] = S[k] * e + 1; }
            p.a[i] = a; p.b[i] = b; T += w;
        }
        p.T = T; return p;
    }
    static void Clip(Cfg c, double[] th)
    {
        th[0] = Math.Max(th[0], c.lo_mu); double sa = 0; for (int k = 0; k < c.K; k++) { th[k + 1] = Math.Max(th[k + 1], c.alMin[k]); sa += th[k + 1] / c.beta[k]; }
        if (sa > c.cap) for (int k = 0; k < c.K; k++) th[k + 1] *= c.cap / sa;
    }
    // quadrature of one gap: dF (at th), dJ (at ths), lg = int lam* ln(lam_th/lam*)
    static void Quad(Cfg c, double[] S, double w, double[] th, bool needF, bool needJ, bool needLg, double[,] dF, double[,] dJ, out double lg)
    {
        int d = c.d, K = c.K; lg = 0; double bmax = c.beta.Max(); int nsub = Math.Max(2, Math.Min(12, (int)Math.Ceiling(w * bmax / 1.5))); double h = w / nsub; var phi = new double[d];
        if (needF) Array.Clear(dF, 0, dF.Length); if (needJ) Array.Clear(dJ, 0, dJ.Length);
        for (int j = 0; j < nsub; j++) for (int q = 0; q < 8; q++)
            {
                double u = h * (j + 0.5 * (gx[q] + 1)), wt = 0.5 * h * gw[q]; phi[0] = 1; for (int k = 0; k < K; k++) phi[k + 1] = S[k] * Math.Exp(-c.beta[k] * u);
                double lt = 0, ls = 0; for (int p = 0; p < d; p++) { lt += phi[p] * th[p]; ls += phi[p] * c.ths[p]; }
                if (needLg) lg += wt * ls * Math.Log(lt / ls);
                for (int p = 0; p < d; p++) for (int r = 0; r <= p; r++) { double vv = phi[p] * phi[r] * wt; if (needF) dF[p, r] += vv / lt; if (needJ) dJ[p, r] += vv / ls; }
            }
    }
    static double Quadf(double[,] A, double[] v, int n) { double s = 0; for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) s += v[i] * A[i, j] * v[j]; return s; }

    // ---- sequential methods; rec(i, th, R) is called after event i (1-based) if i is a checkpoint ----
    static void RunNewton(Cfg c, Path p, bool needR, Action<int, double[], double> rec, int[] cps)
    {
        int d = c.d, K = c.K; var F = (double[,])c.F0.Clone(); var th = (double[])c.th0.Clone(); var dF = new double[d, d]; double R = 0, lg; int gi = 0; var g = new double[d];
        for (int i = 0; i < p.n; i++)
        {
            Quad(c, p.S[i], p.w[i], th, true, false, needR, dF, null, out lg);
            if (needR) { double bt = 0, bs = 0; for (int q = 0; q < d; q++) { bt += p.b[i][q] * th[q]; bs += p.b[i][q] * c.ths[q]; } R += bt - bs - lg; }
            double at = 0; for (int q = 0; q < d; q++) at += p.a[i][q] * th[q];
            for (int q = 0; q < d; q++) for (int r = 0; r <= q; r++) { F[q, r] += dF[q, r]; if (r != q) F[r, q] = F[q, r]; g[q] = p.b[i][q] - p.a[i][q] / at; }
            var st = Solve(F, d, g); for (int q = 0; q < d; q++) th[q] -= st[q]; Clip(c, th);
            if (gi < cps.Length && i + 1 == cps[gi]) { rec(i + 1, th, R); gi++; }
        }
    }
    static void RunHW(Cfg c, Path p, bool needR, char kind, double eta0, Action<int, double[], double> rec, int[] cps)
    {
        int d = c.d; var th = (double[])c.th0.Clone(); double R = 0, lg; int gi = 0; var dummy = new double[d, d];
        for (int i = 0; i < p.n; i++)
        {
            if (needR) { Quad(c, p.S[i], p.w[i], th, false, false, true, dummy, dummy, out lg); double bt = 0, bs = 0; for (int q = 0; q < d; q++) { bt += p.b[i][q] * th[q]; bs += p.b[i][q] * c.ths[q]; } R += bt - bs - lg; }
            double at = 0; for (int q = 0; q < d; q++) at += p.a[i][q] * th[q];
            double eta = kind == 'c' ? eta0 : kind == 's' ? eta0 / Math.Sqrt(i + 1.0) : eta0 / (i + 1.0);
            for (int q = 0; q < d; q++) { double gq = p.b[i][q] - p.a[i][q] / at; th[q] *= Math.Exp(-eta * th[q] * gq); }
            Clip(c, th);
            if (gi < cps.Length && i + 1 == cps[gi]) { rec(i + 1, th, R); gi++; }
        }
    }
    static void RunLS(Cfg c, Path p, bool needR, double rho, Action<int, double[], double> rec, int[] cps)
    {
        int d = c.d, K = c.K; var A = new double[d, d]; var bb = new double[d]; var th = (double[])c.th0.Clone(); double R = 0, lg; int gi = 0; var dummy = new double[d, d];
        var L0 = new double[d, d]; for (int q = 0; q < d; q++) L0[q, q] = c.lam_bar * c.F0[q, q];            // prior in the LS metric (A ~ lam_bar * Fisher)
        var rhs = new double[d]; var M = new double[d, d];
        for (int i = 0; i < p.n; i++)
        {
            if (needR) { Quad(c, p.S[i], p.w[i], th, false, false, true, dummy, dummy, out lg); double bt = 0, bs = 0; for (int q = 0; q < d; q++) { bt += p.b[i][q] * th[q]; bs += p.b[i][q] * c.ths[q]; } R += bt - bs - lg; }
            double w = p.w[i]; var S = p.S[i];
            for (int q = 0; q < d; q++) for (int r = 0; r <= q; r++)
                {
                    A[q, r] *= rho;
                    // int phi_q phi_r du with phi_0 = 1, phi_k = S_k e^{-beta_k u}
                    double val;
                    if (q == 0 && r == 0) val = w;
                    else if (r == 0) val = S[q - 1] * (1 - Math.Exp(-c.beta[q - 1] * w)) / c.beta[q - 1];
                    else val = S[q - 1] * S[r - 1] * (1 - Math.Exp(-(c.beta[q - 1] + c.beta[r - 1]) * w)) / (c.beta[q - 1] + c.beta[r - 1]);
                    A[q, r] += val; if (r != q) A[r, q] = A[q, r];
                }
            for (int q = 0; q < d; q++) bb[q] = rho * bb[q] + p.a[i][q];
            for (int q = 0; q < d; q++) { rhs[q] = bb[q] + L0[q, q] * c.th0[q]; for (int r = 0; r < d; r++) M[q, r] = A[q, r] + (q == r ? L0[q, q] : 0); }
            var sol = Solve(M, d, rhs); for (int q = 0; q < d; q++) th[q] = sol[q]; Clip(c, th);
            if (gi < cps.Length && i + 1 == cps[gi]) { rec(i + 1, th, R); gi++; }
        }
    }
    // ---- batch MLE: Newton with backtracking on all stored gaps up to n ----
    static double Lfun(Path p, int n, double[] t, int d) { double L = 0; for (int j = 0; j < n; j++) { double ap = 0, bp = 0; for (int q = 0; q < d; q++) { ap += p.a[j][q] * t[q]; bp += p.b[j][q] * t[q]; } if (ap <= 1e-12) return double.PositiveInfinity; L += bp - Math.Log(ap); } return L; }
    static double[] BatchFit(Cfg c, Path p, int n, double[] start, out int iters)
    {
        int d = c.d; var tm = (double[])start.Clone(); double Lcur = Lfun(p, n, tm, d); iters = 0;
        for (int it = 0; it < 40; it++)
        {
            var gr = new double[d]; var H = new double[d, d];
            for (int j = 0; j < n; j++) { double ap = 0; for (int q = 0; q < d; q++) ap += p.a[j][q] * tm[q]; double inv = 1 / ap; for (int q = 0; q < d; q++) { gr[q] += p.b[j][q] - p.a[j][q] * inv; for (int r = 0; r <= q; r++) H[q, r] += p.a[j][q] * p.a[j][r] * inv * inv; } }
            for (int q = 0; q < d; q++) { for (int r = 0; r < q; r++) H[r, q] = H[q, r]; H[q, q] += 1e-12; }
            var st = Solve(H, d, gr); double ts = 1; bool ok = false; iters++;
            for (int bt = 0; bt < 30; bt++) { var tn = new double[d]; for (int q = 0; q < d; q++) tn[q] = tm[q] - ts * st[q]; Clip(c, tn); double Ln = Lfun(p, n, tn, d); if (Ln < Lcur - 1e-12 * Math.Abs(Lcur)) { tm = tn; Lcur = Ln; ok = true; break; } ts /= 2; }
            double gn = 0; for (int q = 0; q < d; q++) gn += gr[q] * st[q]; if (!ok || gn < 1e-10) break;
        }
        return tm;
    }

    static Cfg MakeCfg(int K, double[] beta, double mu, double eps)
    {
        int d = K + 1; var c = new Cfg { K = K, d = d, beta = beta, mu_s = mu, eps = eps, rho0 = 1.0 }; c.lam_bar = mu / eps;
        c.al_s = new double[K]; for (int k = 0; k < K; k++) c.al_s[k] = (1 - eps) * beta[k] / K;
        c.ths = new double[d]; c.ths[0] = mu; for (int k = 0; k < K; k++) c.ths[k + 1] = c.al_s[k];
        c.th0 = (double[])c.ths.Clone(); c.th0[0] = 1.3 * mu; for (int k = 0; k < K; k++) c.th0[k + 1] = 0.8 * c.al_s[k];
        c.F0 = new double[d, d]; c.F0[0, 0] = c.rho0 / c.lam_bar; for (int k = 0; k < K; k++) c.F0[k + 1, k + 1] = c.rho0 * c.lam_bar / (beta[k] * beta[k]);
        c.lo_mu = 0.1 * mu; c.cap = 1 - eps / 2; c.alMin = new double[K]; for (int k = 0; k < K; k++) c.alMin[k] = 0.05 * c.al_s[k];
        return c;
    }
    static string F(double x) { return x.ToString("R", CultureInfo.InvariantCulture); }

    static void Main(string[] args)
    {
        var kv = new Dictionary<string, string>(); foreach (var s in args) { var p = s.Split('='); kv[p[0]] = p[1]; }
        Func<string, string, string> get = (k, dflt) => kv.ContainsKey(k) ? kv[k] : dflt; var ci = CultureInfo.InvariantCulture;
        string mode = get("mode", "acc"); int K = int.Parse(get("K", "1")); var beta = get("beta", K == 1 ? "1" : "1,2").Split(',').Select(x => double.Parse(x, ci)).ToArray();
        double mu = double.Parse(get("mu", K == 1 ? "1" : "1.5"), ci), eps = double.Parse(get("eps", "0.1"), ci); int n = int.Parse(get("n", "100000")), reps = int.Parse(get("reps", "32")), seed0 = int.Parse(get("seed", "21000"));
        var c = MakeCfg(K, beta, mu, eps); int d = c.d;
        var cpl = new List<int>(); for (double x = 30; x <= n * 1.0001; x *= Math.Pow(10, 0.25)) cpl.Add((int)Math.Round(x)); if (cpl.Last() != n) cpl.Add(n); var cps = cpl.ToArray();
        var etas = new double[] { 0.003, 0.01, 0.03, 0.1, 0.3, 1.0 };
        if (mode == "acc")
        {
            var outs = new string[reps];
            Parallel.For(0, reps, new ParallelOptions { MaxDegreeOfParallelism = 8 }, r =>
            {
                var path = Gen(c, n, seed0 + r); var sb = new StringBuilder();
                // oracle information at checkpoints
                var Jc = new Dictionary<int, double[,]>(); var Jm = new double[d, d]; var dJ = new double[d, d]; double lg0; int gi = 0; var Tc = new Dictionary<int, double>(); double Tt = 0;
                for (int i = 0; i < n; i++) { Quad(c, path.S[i], path.w[i], c.ths, false, true, false, null, dJ, out lg0); for (int q = 0; q < d; q++) for (int s2 = 0; s2 <= q; s2++) { Jm[q, s2] += dJ[q, s2]; Jm[s2, q] = Jm[q, s2]; } Tt += path.w[i]; if (gi < cps.Length && i + 1 == cps[gi]) { Jc[i + 1] = (double[,])Jm.Clone(); Tc[i + 1] = Tt; gi++; } }
                foreach (var cp in cps) { sb.Append("#I ").Append(r).Append(' ').Append(cp).Append(' ').Append(F(Tc[cp])); for (int q = 0; q < d; q++) for (int s2 = 0; s2 <= q; s2++) sb.Append(' ').Append(F(Jc[cp][q, s2])); sb.Append('\n'); }
                Action<string, int, double[], double> emit = (name, i, th, R) =>
                {
                    var v = new double[d]; for (int q = 0; q < d; q++) v[q] = th[q] - c.ths[q];
                    sb.Append(r).Append(' ').Append(name).Append(' ').Append(i).Append(' ').Append(F(Tc[i])).Append(' ').Append(double.IsNaN(R) ? "nan" : F(R)).Append(' ').Append(F(Quadf(Jc[i], v, d)));
                    for (int q = 0; q < d; q++) sb.Append(' ').Append(F(v[q])); sb.Append('\n');
                };
                RunNewton(c, path, true, (i, th, R) => emit("newton", i, th, R), cps);
                foreach (char kd in new[] { 'c', 's', 'i' }) foreach (double e0 in etas) { char k2 = kd; double ee = e0; RunHW(c, path, true, k2, ee, (i, th, R) => emit("hw_" + k2 + ee.ToString(ci), i, th, R), cps); }
                foreach (double rho in new[] { 1.0, 1 - 1.0 / 2000, 1 - 1.0 / 200 }) { double rr = rho; string nm = rr == 1.0 ? "ls_r1" : rr > 0.999 ? "ls_r2000" : "ls_r200"; RunLS(c, path, true, rr, (i, th, R) => emit(nm, i, th, R), cps); }
                var warm = (double[])c.th0.Clone(); foreach (var cp in cps) { int it; warm = BatchFit(c, path, cp, warm, out it); emit("batch", cp, warm, double.NaN); }
                outs[r] = sb.ToString();
            });
            var w = new System.IO.StreamWriter(Console.OpenStandardOutput());
            w.Write("# K=" + K + " beta=" + string.Join(",", beta) + " mu=" + mu.ToString(ci) + " eps=" + eps.ToString(ci) + " n=" + n + " reps=" + reps + " d=" + d + "\n");
            for (int r = 0; r < reps; r++) w.Write(outs[r]);
            w.Flush();
        }
        else
        {   // timing: single thread, one path, minimum of 3 passes (accuracy not recorded)
            var path = Gen(c, n, seed0); var none = new int[0]; Action<int, double[], double> nop = (i, th, R) => { };
            Console.WriteLine("# timing K=" + K + " eps=" + eps + " n=" + n + " (microseconds per event, minimum of 3 passes; data generation excluded)");
            Func<Action, double> tm = act => { double best = double.MaxValue; for (int t = 0; t < 3; t++) { var sw = Stopwatch.StartNew(); act(); sw.Stop(); best = Math.Min(best, sw.Elapsed.TotalMilliseconds); } return best; };
            double t1 = tm(() => RunNewton(c, path, false, nop, none)); Console.WriteLine(String.Format(ci, "newton        {0,9:F3} us/event   ({1:F1} ms total)", 1000 * t1 / n, t1));
            double t2 = tm(() => RunHW(c, path, false, 's', 0.1, nop, none)); Console.WriteLine(String.Format(ci, "hw (mirror)   {0,9:F3} us/event   ({1:F1} ms total)", 1000 * t2 / n, t2));
            double t3 = tm(() => RunLS(c, path, false, 1.0, nop, none)); Console.WriteLine(String.Format(ci, "ls (EWLS)     {0,9:F3} us/event   ({1:F1} ms total)", 1000 * t3 / n, t3));
            Console.WriteLine("batch MLE refit time vs n (cold start from the near start, then one warm-started refit):");
            foreach (int m in new[] { 1000, 10000, 100000, Math.Min(n, 200000) }.Distinct().Where(x => x <= n))
            {
                int it = 0; double tc = tm(() => { BatchFit(c, path, m, c.th0, out it); }); var sol = BatchFit(c, path, m, c.th0, out it);
                int it2 = 0; double tw = tm(() => { BatchFit(c, path, m, sol, out it2); });
                Console.WriteLine(String.Format(ci, "  n={0,7}: cold refit {1,9:F2} ms ({2} Newton its)   warm refit (already at optimum) {3,8:F2} ms   -> {4:F3} us per event of data per refit", m, tc, it, tw, 1000 * tc / m));
            }
            var warmt = (double[])c.th0.Clone(); var sw2 = Stopwatch.StartNew(); foreach (var cp in cps) { int it; warmt = BatchFit(c, path, cp, warmt, out it); } sw2.Stop();
            Console.WriteLine(String.Format(ci, "batch refits at the {0} geometric checkpoints (warm-started chain): total {1:F1} ms = {2:F3} us/event amortised over n = {3}", cps.Length, sw2.Elapsed.TotalMilliseconds, 1000 * sw2.Elapsed.TotalMilliseconds / n, n));
            int itn; double tn = tm(() => { int it; BatchFit(c, path, n, c.th0, out it); }); BatchFit(c, path, n, c.th0, out itn);
            Console.WriteLine(String.Format(ci, "keeping the batch estimate current after EVERY event with one Newton step on all past data: cost per event = (time of one Newton iteration at size k) ~ {0:F1} us at k = n = {1}, average over k <= n: {2:F1} us/event", 1000 * tn / itn, n, 1000 * tn / itn / 2));
        }
    }
}
