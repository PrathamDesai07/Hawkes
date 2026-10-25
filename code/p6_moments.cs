// Problem 6 (D): stationary negative moments E_pi[lam^-q] of the K-exponential Hawkes intensity (time average over a long run, GL quadrature inside each gap),
// against the rigorous recursion bound E w^-(p+1) <= eps*r/(mu-(p+1)*B*n^2/2) * E w^-p (w = mu + b Z, Z = sum alpha_k S_k/beta_k, lam >= w, r = B/b).
// Also the stationary information I = E[a a'/lam] (time average of phi phi'/lam) for det I / eigenvalue checks.
// Build: csc /o+ /out:p6_moments.exe p6_moments.cs   Run: p6_moments.exe K=1 beta=1 mu=1 eps=0.1 tmax=200000 reps=8
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Threading.Tasks;

static class P6M
{
    static readonly double[] gx = { -0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498, 0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363 };
    static readonly double[] gw = { 0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620, 0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763 };
    static readonly double[] qs = { 1.0, 1.5, 2.0, 2.5, 3.0 };

    static double[] Run(int K, double[] beta, double mu, double eps, double tmax, int seed, double warmT, out double[,] I)
    {
        int d = K + 1; var rng = new Random(seed);
        var al = new double[K]; for (int k = 0; k < K; k++) al[k] = (1 - eps) * beta[k] / K;
        var S = new double[K]; double Tt = 0, T0 = -1; var acc = new double[qs.Length]; I = new double[d, d];
        double bmax = beta.Max(); var phi = new double[d];
        while (Tt < warmT + tmax)
        {
            double lam0 = mu; for (int k = 0; k < K; k++) lam0 += al[k] * S[k];
            double w = 0;
            while (true)
            {
                double u = rng.NextDouble(); while (u <= 0) u = rng.NextDouble();
                w += -Math.Log(u) / lam0;
                double l = mu; for (int k = 0; k < K; k++) l += al[k] * S[k] * Math.Exp(-beta[k] * w);
                if (rng.NextDouble() * lam0 <= l) break;
            }
            if (Tt >= warmT)
            {
                int nsub = Math.Max(1, Math.Min(12, (int)Math.Ceiling(w * bmax / 1.5))); double h = w / nsub;
                for (int j = 0; j < nsub; j++) for (int q = 0; q < 8; q++)
                    {
                        double uu = h * (j + 0.5 * (gx[q] + 1)), wt = 0.5 * h * gw[q];
                        phi[0] = 1; double lam = mu; for (int k = 0; k < K; k++) { phi[k + 1] = S[k] * Math.Exp(-beta[k] * uu); lam += al[k] * phi[k + 1]; }
                        for (int a = 0; a < qs.Length; a++) acc[a] += wt * Math.Pow(lam, -qs[a]);
                        for (int p = 0; p < d; p++) for (int r2 = 0; r2 <= p; r2++) I[p, r2] += wt * phi[p] * phi[r2] / lam;
                    }
            }
            for (int k = 0; k < K; k++) S[k] = S[k] * Math.Exp(-beta[k] * w) + 1;
            Tt += w;
        }
        double T = Tt - warmT;
        for (int a = 0; a < acc.Length; a++) acc[a] /= T;
        for (int p = 0; p < d; p++) for (int r2 = 0; r2 <= p; r2++) { I[p, r2] /= T; I[r2, p] = I[p, r2]; }
        return acc;
    }

    static double[] Eig(double[,] A, int d)
    {
        var M = (double[,])A.Clone();
        for (int sweep = 0; sweep < 60; sweep++)
        {
            double off = 0; for (int i = 0; i < d; i++) for (int j = 0; j < i; j++) off += M[i, j] * M[i, j];
            if (off < 1e-26) break;
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
        double mu = double.Parse(get("mu", "1"), ci), eps = double.Parse(get("eps", "0.1"), ci), tmax = double.Parse(get("tmax", "100000"), ci);
        int reps = int.Parse(get("reps", "8")); int d = K + 1;
        var accs = new double[reps][]; var Is = new double[reps][,];
        Parallel.For(0, reps, new ParallelOptions { MaxDegreeOfParallelism = 8 }, r => { double[,] I; accs[r] = Run(K, beta, mu, eps, tmax, 900 + r, 20.0 / (beta[0] * eps), out I); Is[r] = I; });
        double b = beta.Min(), B = beta.Max(), n = 1 - eps, rr = B / b, kEff = 2 * mu / (B * n * n);
        Console.WriteLine("# K=" + K + " beta=" + string.Join(",", beta) + " mu=" + mu + " eps=" + eps + " k_eff=2mu/(B n^2)=" + kEff.ToString("F3", ci));
        // recursion bound: E w^-(p+1) <= c_p E w^-p, c_p = eps r/(mu-(p+1)B n^2/2); fractional: E w^-p <= (E w^-1)^p for p in [0,1]
        Func<double, double> c = p => (mu - (p + 1) * B * n * n / 2) > 0 ? eps * rr / (mu - (p + 1) * B * n * n / 2) : double.NaN;
        Func<double, double> bound = q =>
        {
            double v = 1; double qq = q;
            while (qq > 1 + 1e-12) { double p = qq - 1; v *= c(p); qq -= 1; if (double.IsNaN(v)) return double.NaN; }
            // now qq in (0,1]: E w^-qq <= (E w^-1)^qq, with E w^-1 <= c_0
            double c0 = c(0); if (double.IsNaN(c0)) return double.NaN;
            double rest = Math.Pow(c0, qq);
            // last step for fractional qq<1 uses E w^-(p+1) with p = qq-1 <0: just Jensen
            return v * rest;
        };
        Console.WriteLine("#  q   MC E[lam^-q] (se)     bound   bound/MC    lam_bar^-q*C?  (MC*(lam_bar)^q)");
        double lamBar = mu / eps;
        for (int a = 0; a < qs.Length; a++)
        {
            double m = accs.Average(x => x[a]); double se = Math.Sqrt(accs.Sum(x => (x[a] - m) * (x[a] - m)) / Math.Max(1, reps - 1) / reps);
            double bd = bound(qs[a]);
            Console.WriteLine(String.Format(ci, "{0,5:F1} {1,12:E4} ({2,9:E2}) {3,12:E4} {4,8:F2} {5,12:F4}", qs[a], m, se, bd, bd / m, m * Math.Pow(lamBar, qs[a])));
        }
        var Ib = new double[d, d]; for (int p = 0; p < d; p++) for (int q = 0; q <= p; q++) { double v = Is.Average(x => x[p, q]); Ib[p, q] = v; Ib[q, p] = v; }
        var ev = Eig(Ib, d); double det = 1; foreach (var e in ev) det *= e;
        Console.WriteLine("# I eigenvalues: " + string.Join(" ", ev.Select(e => e.ToString("E4", ci))) + "   det I = " + det.ToString("E4", ci) + "   eps*lam_max = " + (ev[d - 1] * eps).ToString("F3", ci) + "   lam_min/eps = " + (ev[0] / eps).ToString("F4", ci));
    }
}
