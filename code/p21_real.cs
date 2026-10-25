// Problem 21, real-data half: baselines on the Deribit BTC option-trade stream (code/p20_data/events.bin; windows 13-20 UTC, 365 days; K = 3, beta = (0.1, 1, 10)/s).
// Protocol: every method predicts BEFORE it sees each event (prequential), so the per-event point-process log-loss l_i(theta_i) = b_i'theta_i - ln a_i'theta_i (+ the final partial gap b_end'theta)
// is an out-of-sample score.  Methods run from day 0; day-level losses are written so that any evaluation period can be chosen afterwards.
//   newton_m5 / newton_m1 / newton_g : our online Newton, forgetting eta = 1/(m N_prevday) (m = 5, 1) or growing window (eta = 0)
//   hw_{c|s|i}<eta0> : Hall-Willett-type mirror descent on the log-parameters (constant / eta0/sqrt(i) / eta0/i)
//   ls_r1 / ls_r20000 / ls_r2000 / ls_r200 : Solo-type exponentially weighted least squares (rho per event)
//   batch_all : theta refit (batch MLE, warm start) on ALL previous days at the start of each day, held fixed during the day
//   batch_30d : same on the previous 30 days only (a standard remedy for non-stationarity)
//   oracle_day: in-sample MLE of each day (best per-day parameter; a tracking lower bound, not a predictor)
//   pooled_all / pooled_test : best fixed parameter in hindsight over all days / over days >= test0 (in-sample)
// Output lines:  LOSS day nev method loss        DIST day method e_dist dn_hat      TIME method us_per_event
// Build: csc /o+ /out:p21_real.exe p21_real.cs     Run: p21_real.exe mode=run   |   mode=time
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;

static class P21R
{
    static double[] ev; static CultureInfo ci = CultureInfo.InvariantCulture;
    static double[] beta = { 0.1, 1, 10 }; static int K = 3, d = 4;
    static int LowerBound(double t) { int lo = 0, hi = ev.Length; while (lo < hi) { int m = (lo + hi) / 2; if (ev[m] < t) lo = m + 1; else hi = m; } return lo; }
    class Win { public int n; public double t0, t1; public double[] w; public double[][] a, b, S0; public double[] bEnd; }
    static Win Build(double t0, double t1)
    {
        int i0 = LowerBound(t0), i1 = LowerBound(t1), iw = LowerBound(t0 - 600); var S = new double[K]; double tp = t0 - 600;
        for (int i = iw; i < i0; i++) { double w = ev[i] - tp; for (int k = 0; k < K; k++) S[k] = S[k] * Math.Exp(-beta[k] * w) + 1; tp = ev[i]; }
        { double w = t0 - tp; for (int k = 0; k < K; k++) S[k] *= Math.Exp(-beta[k] * w); }
        var W = new Win { n = i1 - i0, t0 = t0, t1 = t1, w = new double[i1 - i0], a = new double[i1 - i0][], b = new double[i1 - i0][], S0 = new double[i1 - i0][] }; double last = t0;
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
    static void Clip(double[] th) { th[0] = Math.Max(th[0], 1e-6); double sa = 0; for (int k = 0; k < K; k++) { th[k + 1] = Math.Max(th[k + 1], 1e-8); sa += th[k + 1] / beta[k]; } if (sa > 0.999) for (int k = 0; k < K; k++) th[k + 1] *= 0.999 / sa; }
    static double Nhat(double[] th) { double n = 0; for (int k = 0; k < K; k++) n += th[k + 1] / beta[k]; return n; }
    static double Nll(List<Win> ws, double[] th)
    {
        double L = 0; foreach (var W in ws) { for (int j = 0; j < W.n; j++) { double ap = 0, bp = 0; for (int q = 0; q < d; q++) { ap += W.a[j][q] * th[q]; bp += W.b[j][q] * th[q]; } if (ap <= 1e-300) return double.PositiveInfinity; L += bp - Math.Log(ap); } for (int q = 0; q < d; q++) L += W.bEnd[q] * th[q]; }
        return L;
    }
    static double[] Fit(List<Win> ws, double[] start, out double[,] H, out double nll)
    {
        var tm = (double[])start.Clone(); Clip(tm); nll = Nll(ws, tm); H = new double[d, d];
        for (int it = 0; it < 80; it++)
        {
            var gr = new double[d]; H = new double[d, d];
            foreach (var W in ws) { for (int j = 0; j < W.n; j++) { double ap = 0; for (int q = 0; q < d; q++) ap += W.a[j][q] * tm[q]; double inv = 1 / ap; for (int q = 0; q < d; q++) { gr[q] += W.b[j][q] - W.a[j][q] * inv; for (int r = 0; r <= q; r++) H[q, r] += W.a[j][q] * W.a[j][r] * inv * inv; } } for (int q = 0; q < d; q++) gr[q] += W.bEnd[q]; }
            for (int q = 0; q < d; q++) { for (int r = 0; r < q; r++) H[r, q] = H[q, r]; H[q, q] += 1e-9 * (1 + H[q, q]); }
            var st = Solve(H, d, gr); double ts = 1; bool ok = false;
            for (int bt = 0; bt < 40; bt++) { var tn = new double[d]; for (int q = 0; q < d; q++) tn[q] = tm[q] - ts * st[q]; Clip(tn); double Ln = Nll(ws, tn); if (Ln < nll - 1e-10 * Math.Abs(nll)) { tm = tn; nll = Ln; ok = true; break; } ts /= 2; }
            if (!ok) break;
        }
        return tm;
    }
    static readonly double[] gx = { -0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498, 0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363 };
    static readonly double[] gw = { 0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620, 0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763 };

    abstract class Est
    {
        public string name; public double[] th;
        public virtual void StartDay(Win W, int day, double nPrev) { }
        public abstract void Update(Win W, int j, long gi);     // uses the event after the prediction; th is the prediction for event j before the call
    }
    class Newton : Est
    {
        double m; double[,] F; double eta; double[,] dF = new double[4, 4]; double[] g = new double[4]; bool init = false; double rho0 = 1;
        public Newton(string nm, double mult, double[] th0) { name = nm; m = mult; th = (double[])th0.Clone(); }
        public override void StartDay(Win W, int day, double nPrev)
        {
            if (!init) { double lb = W.n / (W.t1 - W.t0); F = new double[d, d]; F[0, 0] = rho0 / lb; for (int k = 0; k < K; k++) F[k + 1, k + 1] = rho0 * lb / (beta[k] * beta[k]); init = true; }
            eta = m > 0 ? 1.0 / (m * nPrev) : 0.0;
        }
        public override void Update(Win W, int j, long gi)
        {
            Array.Clear(dF, 0, dF.Length); double w = W.w[j]; int nsub = Math.Max(2, Math.Min(12, (int)Math.Ceiling(w * beta.Max() / 1.5))); double h = w / nsub; var phi = new double[d];
            for (int s = 0; s < nsub; s++) for (int q = 0; q < 8; q++) { double u = h * (s + 0.5 * (gx[q] + 1)), wt = 0.5 * h * gw[q]; phi[0] = 1; for (int k = 0; k < K; k++) phi[k + 1] = W.S0[j][k] * Math.Exp(-beta[k] * u); double lt = 0; for (int p = 0; p < d; p++) lt += phi[p] * th[p]; for (int p = 0; p < d; p++) for (int r = 0; r <= p; r++) dF[p, r] += wt * phi[p] * phi[r] / lt; }
            double at = 0; for (int p = 0; p < d; p++) at += W.a[j][p] * th[p];
            for (int p = 0; p < d; p++) { for (int r = 0; r <= p; r++) { F[p, r] = (1 - eta) * F[p, r] + dF[p, r]; if (r != p) F[r, p] = F[p, r]; } g[p] = W.b[j][p] - W.a[j][p] / at; }
            var st = Solve(F, d, g); for (int p = 0; p < d; p++) th[p] -= st[p]; Clip(th);
        }
    }
    class HW : Est
    {
        char kind; double eta0;
        public HW(string nm, char k, double e0, double[] th0) { name = nm; kind = k; eta0 = e0; th = (double[])th0.Clone(); }
        public override void Update(Win W, int j, long gi)
        {
            double at = 0; for (int p = 0; p < d; p++) at += W.a[j][p] * th[p]; double eta = kind == 'c' ? eta0 : kind == 's' ? eta0 / Math.Sqrt(gi + 1.0) : eta0 / (gi + 1.0);
            for (int p = 0; p < d; p++) { double gq = W.b[j][p] - W.a[j][p] / at; th[p] *= Math.Exp(-eta * th[p] * gq); } Clip(th);
        }
    }
    class LS : Est
    {
        double rho; double[,] A = new double[4, 4]; double[] bb = new double[4]; double[,] L0 = new double[4, 4]; double[] th0; bool init = false; double[,] M = new double[4, 4]; double[] rhs = new double[4];
        public LS(string nm, double r, double[] t0) { name = nm; rho = r; th0 = (double[])t0.Clone(); th = (double[])t0.Clone(); }
        public override void StartDay(Win W, int day, double nPrev) { if (!init) { double lb = W.n / (W.t1 - W.t0); for (int q = 0; q < d; q++) L0[q, q] = lb * (q == 0 ? 1 / lb : lb / (beta[q - 1] * beta[q - 1])); init = true; } }
        public override void Update(Win W, int j, long gi)
        {
            double w = W.w[j]; var S = W.S0[j];
            for (int q = 0; q < d; q++) for (int r = 0; r <= q; r++)
                {
                    double val; if (q == 0 && r == 0) val = w; else if (r == 0) val = S[q - 1] * (1 - Math.Exp(-beta[q - 1] * w)) / beta[q - 1]; else val = S[q - 1] * S[r - 1] * (1 - Math.Exp(-(beta[q - 1] + beta[r - 1]) * w)) / (beta[q - 1] + beta[r - 1]);
                    A[q, r] = rho * A[q, r] + val; if (r != q) A[r, q] = A[q, r];
                }
            for (int q = 0; q < d; q++) { bb[q] = rho * bb[q] + W.a[j][q]; rhs[q] = bb[q] + L0[q, q] * th0[q]; for (int r = 0; r < d; r++) M[q, r] = A[q, r] + (q == r ? L0[q, q] : 0); }
            var sol = Solve(M, d, rhs); for (int q = 0; q < d; q++) th[q] = sol[q]; Clip(th);
        }
    }
    static double DayLoss(Win W, double[] th)
    {
        double L = 0; for (int j = 0; j < W.n; j++) { double ap = 0, bp = 0; for (int q = 0; q < d; q++) { ap += W.a[j][q] * th[q]; bp += W.b[j][q] * th[q]; } L += bp - Math.Log(ap); } for (int q = 0; q < d; q++) L += W.bEnd[q] * th[q]; return L;
    }

    static void Main(string[] args)
    {
        var kv = new Dictionary<string, string>(); foreach (var s in args) { int p = s.IndexOf('='); kv[s.Substring(0, p)] = s.Substring(p + 1); }
        Func<string, string, string> get = (k, dflt) => kv.ContainsKey(k) ? kv[k] : dflt;
        var bytes = File.ReadAllBytes(get("data", @"p20_data\events.bin")); ev = new double[bytes.Length / 8]; Buffer.BlockCopy(bytes, 0, ev, 0, bytes.Length);
        string mode = get("mode", "run"); int nd = int.Parse(get("days", "365")); int test0 = int.Parse(get("test0", "60")); double h0 = 13, h1 = 20;
        var wins = new List<Win>(); for (int dd = 0; dd < nd; dd++) wins.Add(Build(dd * 86400.0 + h0 * 3600, dd * 86400.0 + h1 * 3600));
        var th0 = new double[d]; { double lb = wins[0].n / (wins[0].t1 - wins[0].t0); th0[0] = 0.5 * lb; for (int k = 0; k < K; k++) th0[k + 1] = 0.3 * beta[k] / K; }
        var etas = new double[] { 0.003, 0.01, 0.03, 0.1, 0.3, 1.0 };
        Func<List<Est>> make = () =>
        {
            var l = new List<Est>(); l.Add(new Newton("newton_m5", 5, th0)); l.Add(new Newton("newton_m1", 1, th0)); l.Add(new Newton("newton_g", 0, th0)); if (get("more", "0") == "1" || get("more", "0") == "2") foreach (double mm in new[] { 20.0, 2.0, 0.5, 0.2, 0.1 }) l.Add(new Newton("newton_m" + mm.ToString(ci), mm, th0)); if (get("more", "0") == "2") { l.Clear(); foreach (double mm in new[] { 0.05, 0.02, 0.01, 0.005 }) l.Add(new Newton("newton_m" + mm.ToString(ci), mm, th0)); foreach (char kd in new[] { 'c', 's', 'i' }) foreach (double e0 in new[] { 3.0, 10.0 }) l.Add(new HW("hw_" + kd + e0.ToString(ci), kd, e0, th0)); foreach (int q in new[] { 100, 50, 20 }) l.Add(new LS("ls_r" + q, 1 - 1.0 / q, th0)); return l; }
            foreach (char kd in new[] { 'c', 's', 'i' }) foreach (double e0 in etas) l.Add(new HW("hw_" + kd + e0.ToString(ci), kd, e0, th0));
            l.Add(new LS("ls_r1", 1.0, th0)); l.Add(new LS("ls_r20000", 1 - 1.0 / 20000, th0)); l.Add(new LS("ls_r2000", 1 - 1.0 / 2000, th0)); l.Add(new LS("ls_r200", 1 - 1.0 / 200, th0)); return l;
        };
        if (mode == "run")
        {
            var ests = make(); var thB = (double[])th0.Clone(); var thB30 = (double[])th0.Clone(); long gi = 0; double nPrev = wins[0].n; var distDays = new HashSet<int> { 5, 10, 20, 30, 60, 90, 120, 180, 240, 300, 364 };
            var lossTab = new Dictionary<string, double>();
            for (int dd = 0; dd < nd; dd++)
            {
                var W = wins[dd]; if (W.n < 1) continue;
                // batch predictors: refit on previous days (warm start), fixed during the day
                if (dd > 0 && get("fast", "0") != "1") { double[,] H; double nl; thB = Fit(wins.GetRange(0, dd), thB, out H, out nl); int s0 = Math.Max(0, dd - 30); thB30 = Fit(wins.GetRange(s0, dd - s0), thB30, out H, out nl); }
                if (get("fast", "0") != "1") { Console.WriteLine(String.Format(ci, "LOSS {0} {1} batch_all {2:R}", dd, W.n, DayLoss(W, thB))); Console.WriteLine(String.Format(ci, "LOSS {0} {1} batch_30d {2:R}", dd, W.n, DayLoss(W, thB30))); }
                { double[,] H; double nl; var one = new List<Win> { W }; var tho = Fit(one, thB, out H, out nl); Console.WriteLine(String.Format(ci, "LOSS {0} {1} oracle_day {2:R}", dd, W.n, nl)); }
                foreach (var e in ests)
                {
                    e.StartDay(W, dd, nPrev); double L = 0;
                    for (int j = 0; j < W.n; j++) { double ap = 0, bp = 0; for (int q = 0; q < d; q++) { ap += W.a[j][q] * e.th[q]; bp += W.b[j][q] * e.th[q]; } L += bp - Math.Log(ap); e.Update(W, j, gi + j); }
                    for (int q = 0; q < d; q++) L += W.bEnd[q] * e.th[q];
                    Console.WriteLine(String.Format(ci, "LOSS {0} {1} {2} {3:R}", dd, W.n, e.name, L));
                }
                gi += W.n; nPrev = W.n;
                if (distDays.Contains(dd) && get("fast", "0") != "1")
                {
                    double[,] H; double nl; var thP = Fit(wins.GetRange(0, dd + 1), thB, out H, out nl);
                    foreach (var e in ests) { var df = new double[d]; for (int q = 0; q < d; q++) df[q] = e.th[q] - thP[q]; double ed = 0; for (int p = 0; p < d; p++) for (int q = 0; q < d; q++) ed += df[p] * H[p, q] * df[q]; Console.WriteLine(String.Format(ci, "DIST {0} {1} {2:R} {3:R}", dd, e.name, ed, Nhat(e.th) - Nhat(thP))); }
                    Console.WriteLine(String.Format(ci, "DIST {0} batch_pooled 0 {1:R}", dd, Nhat(thP)));
                }
            }
            // in-sample best fixed parameters in hindsight
            double[,] H2; double nl2; var thAll = Fit(wins, thB, out H2, out nl2); var wt = wins.GetRange(test0, nd - test0); var thT = Fit(wt, thAll, out H2, out nl2);
            for (int dd = 0; dd < nd; dd++) { Console.WriteLine(String.Format(ci, "LOSS {0} {1} pooled_all {2:R}", dd, wins[dd].n, DayLoss(wins[dd], thAll))); Console.WriteLine(String.Format(ci, "LOSS {0} {1} pooled_test {2:R}", dd, wins[dd].n, DayLoss(wins[dd], thT))); }
            Console.WriteLine("THETA pooled_all " + string.Join(" ", thAll.Select(x => x.ToString("G6", ci))) + " nhat " + Nhat(thAll).ToString("F4", ci)); Console.WriteLine("THETA pooled_test " + string.Join(" ", thT.Select(x => x.ToString("G6", ci))) + " nhat " + Nhat(thT).ToString("F4", ci));
        }
        else
        {   // timing: one pass over all windows per method, no loss bookkeeping
            long total = wins.Sum(w => (long)w.n); Console.WriteLine("# timing on the real stream: " + total + " events, " + nd + " windows (microseconds per event; one thread)");
            foreach (var proto in make().Where(e => e.name == "newton_m5" || e.name == "hw_s0.1" || e.name == "ls_r1" || e.name == "ls_r2000"))
            {
                var e = make().First(x => x.name == proto.name); var sw = Stopwatch.StartNew(); long gi = 0; double nPrev = wins[0].n;
                for (int dd = 0; dd < nd; dd++) { var W = wins[dd]; e.StartDay(W, dd, nPrev); for (int j = 0; j < W.n; j++) e.Update(W, j, gi + j); gi += W.n; nPrev = W.n; }
                sw.Stop(); Console.WriteLine(String.Format(ci, "TIME {0} {1:F3} us/event ({2:F0} ms)", e.name, 1000 * sw.Elapsed.TotalMilliseconds / total, sw.Elapsed.TotalMilliseconds));
            }
            var swb = Stopwatch.StartNew(); var thb = (double[])th0.Clone(); for (int dd = 1; dd < nd; dd++) { double[,] H; double nl; thb = Fit(wins.GetRange(0, dd), thb, out H, out nl); } swb.Stop();
            Console.WriteLine(String.Format(ci, "TIME batch_all (daily warm refit on all previous days) {0:F3} us/event ({1:F0} ms total for {2} refits; events seen by refits grow linearly)", 1000 * swb.Elapsed.TotalMilliseconds / total, swb.Elapsed.TotalMilliseconds, nd - 1));
            var swc = Stopwatch.StartNew(); { double[,] H; double nl; Fit(wins, th0, out H, out nl); } swc.Stop(); Console.WriteLine(String.Format(ci, "TIME batch cold refit on the full year: {0:F0} ms = {1:F3} us per data event", swc.Elapsed.TotalMilliseconds, 1000 * swc.Elapsed.TotalMilliseconds / total));
        }
    }
}
