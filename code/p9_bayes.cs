// Problem 9, K=1: Bayes regret of the online Newton predictor over the stability band, the matrix-covariance lower bound of Prop 9.1 evaluated with the
// MEASURED covariance, the Gaussian-experiment benchmark I_G(T), and the exponential-dictionary two-point bound (Prop 9.0).
//   mode=bayes : theta* ~ Unif( mu in mu_c(1 +- w) , alpha in beta(1 - 3eps/2 , 1 - eps/2) ),  S_0 stationary under theta* (observed),
//                F_0 = rho0 * diag(12/width_mu^2, 12/width_alpha^2) (prior precision of the uniform), start at the box centre.
//                at each checkpoint T: E R_T (Newton), s.e., Cov = E[(theta_hat - theta*)(.)'] (2x2), mean J_T/T (2x2), L_emp = ln|Theta_0| - 0.5 ln((2 pi e)^2 det Cov),
//                I_G(T) = 0.5 ln det(I + T S^{1/2} Ibar S^{1/2}),  S = diag(width^2/12).
//   mode=twopoint: fixed theta, constant predictor theta' = theta + delta*u: E R_T / (T u'Iu delta^2 / 2), with ratio range for the c~ bounds.
// Build: csc /o+ /out:p9_bayes.exe p9_bayes.cs     Run: p9_bayes.exe mode=bayes eps=0.1 mu=1 beta=1 reps=400 tmax=3000
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Threading.Tasks;

static class P9
{
    static readonly double[] gx = { -0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498, 0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363 };
    static readonly double[] gw = { 0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620, 0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763 };
    static double U(Random r) { double u; do { u = r.NextDouble(); } while (u <= 0); return u; }
    // closed-form Fisher increments for K=1: (int 1/lam, int s e^{-bu}/lam, int (s e^{-bu})^2/lam)
    static void FI(double mu, double al, double beta, double s, double w, out double i0, out double i1, out double i2)
    {
        double c = al * s, e = Math.Exp(-beta * w);
        i0 = (w + Math.Log((mu + c * e) / (mu + c)) / beta) / mu; double Lam = mu * w + c * (1 - e) / beta;
        i1 = (w - mu * i0) / al; i2 = (Lam - 2 * mu * w + mu * mu * i0) / (al * al);
    }
    static double Gap(Random rng, double mu_s, double al_s, double beta, double s)
    {
        double lam0 = mu_s + al_s * s, w = 0;
        while (true) { w += -Math.Log(U(rng)) / lam0; if (rng.NextDouble() * lam0 <= mu_s + al_s * s * Math.Exp(-beta * w)) return w; }
    }
    class R1 { public double[] R, Tt, vmu, val, J0, J1, J2, n; }

    static R1 Bayes(int seed, double mu_c, double w, double beta, double eps, double rho0, double[] grid, int warm)
    {
        var rng = new Random(seed);
        double wm = 2 * w * mu_c, wa = eps * beta;
        double mu_s = mu_c * (1 - w + 2 * w * rng.NextDouble()), al_s = beta * (1 - 1.5 * eps + eps * rng.NextDouble());
        double mu = mu_c, al = beta * (1 - eps);                     // start at the box centre
        double F00 = rho0 * 12 / (wm * wm), F11 = rho0 * 12 / (wa * wa), F01 = 0;
        double s = 0; for (int i = 0; i < warm; i++) { double ww = Gap(rng, mu_s, al_s, beta, s); s = s * Math.Exp(-beta * ww) + 1; }
        int G = grid.Length, gi = 0; var o = new R1 { R = new double[G], Tt = new double[G], vmu = new double[G], val = new double[G], J0 = new double[G], J1 = new double[G], J2 = new double[G], n = new double[G] };
        double R = 0, Tt = 0, j0 = 0, j1 = 0, j2 = 0; int cnt = 0;
        while (gi < G)
        {
            double wgap = Gap(rng, mu_s, al_s, beta, s); cnt++;
            double e = Math.Exp(-beta * wgap), a1 = s * e, b1 = s * (1 - e) / beta;
            double Lt = mu * wgap + al * b1, Ls = mu_s * wgap + al_s * b1, lg = 0;
            for (int q = 0; q < 8; q++) { double u = 0.5 * wgap * (gx[q] + 1), su = s * Math.Exp(-beta * u); double ls = mu_s + al_s * su, lt = mu + al * su; lg += 0.5 * wgap * gw[q] * ls * Math.Log(lt / ls); }
            R += Lt - Ls - lg;
            double i0, i1, i2; FI(mu_s, al_s, beta, s, wgap, out i0, out i1, out i2); j0 += i0; j1 += i1; j2 += i2;
            double f0, f1, f2; FI(mu, al, beta, s, wgap, out f0, out f1, out f2); F00 += f0; F01 += f1; F11 += f2;
            double g0 = wgap - 1.0 / (mu + al * a1), g1 = b1 - a1 / (mu + al * a1);
            double det = F00 * F11 - F01 * F01; double d0 = (F11 * g0 - F01 * g1) / det, d1 = (-F01 * g0 + F00 * g1) / det;
            mu = Math.Max(mu - d0, 0.1 * mu_c); al = Math.Min(Math.Max(al - d1, 0.02 * beta), beta * (1 - eps / 4));   // clip: stability band relaxed to eps/4
            s = a1 + 1; Tt += wgap;
            while (gi < G && Tt >= grid[gi]) { o.R[gi] = R; o.Tt[gi] = Tt; o.vmu[gi] = mu - mu_s; o.val[gi] = al - al_s; o.J0[gi] = j0; o.J1[gi] = j1; o.J2[gi] = j2; o.n[gi] = cnt; gi++; }
        }
        return o;
    }

    static void Main(string[] args)
    {
        var kv = new Dictionary<string, string>(); foreach (var a in args) { var p = a.Split('='); kv[p[0]] = p[1]; }
        Func<string, string, string> get = (k, d) => kv.ContainsKey(k) ? kv[k] : d; var ci = CultureInfo.InvariantCulture;
        string mode = get("mode", "bayes"); double eps = double.Parse(get("eps", "0.1"), ci), mu = double.Parse(get("mu", "1"), ci), beta = double.Parse(get("beta", "1"), ci);
        double w = double.Parse(get("w", "0.25"), ci), rho0 = double.Parse(get("rho0", "1"), ci), tmax = double.Parse(get("tmax", "3000"), ci), tmin = double.Parse(get("tmin", "10"), ci);
        int reps = int.Parse(get("reps", "200")), seed0 = int.Parse(get("seed", "9000"));
        var grid = new List<double>(); for (double t = tmin; t <= tmax * 1.0001; t *= Math.Pow(10, 0.25)) grid.Add(t); var G = grid.ToArray();
        if (mode == "bayes")
        {
            int warm = (int)(20.0 / (eps * eps)); var all = new R1[reps];
            Parallel.For(0, reps, new ParallelOptions { MaxDegreeOfParallelism = 8 }, r => { all[r] = Bayes(seed0 + r, mu, w, beta, eps, rho0, G, warm); });
            double wm = 2 * w * mu, wa = eps * beta, lnVol = Math.Log(wm * wa);
            Console.WriteLine("# bayes eps=" + eps + " mu=" + mu + " beta=" + beta + " w=" + w + " rho0=" + rho0 + " reps=" + reps + " ln|Theta0|=" + lnVol.ToString("F4", ci) + "  (width_mu=" + wm.ToString("F3", ci) + " width_alpha=" + wa.ToString("F4", ci) + ")");
            Console.WriteLine("# T  E_R  se  L_emp  gap=E_R-L_emp  I_G  E_R-I_G  Cov_mm Cov_aa Cov_ma  kappaeff=T*lmax(I^1/2 Cov I^1/2)/1  lnDetI  lnDetCov*T^2  n_events");
            for (int g = 0; g < G.Length; g++)
            {
                double mR = all.Average(x => x.R[g]), se = Math.Sqrt(all.Sum(x => (x.R[g] - mR) * (x.R[g] - mR)) / Math.Max(1, reps - 1) / reps), mT = all.Average(x => x.Tt[g]);
                double c00 = all.Average(x => x.vmu[g] * x.vmu[g]), c11 = all.Average(x => x.val[g] * x.val[g]), c01 = all.Average(x => x.vmu[g] * x.val[g]);
                double detC = c00 * c11 - c01 * c01;
                double L = lnVol - 0.5 * (2 * Math.Log(2 * Math.PI * Math.E) + Math.Log(detC));
                double I00 = all.Average(x => x.J0[g]) / mT, I01 = all.Average(x => x.J1[g]) / mT, I11 = all.Average(x => x.J2[g]) / mT;
                // I_G = 0.5 ln det(I + T S^{1/2} I S^{1/2}),  S = diag(wm^2/12, wa^2/12)
                double s0 = wm * wm / 12, s1 = wa * wa / 12; double m00 = 1 + mT * s0 * I00, m11 = 1 + mT * s1 * I11, m01 = mT * Math.Sqrt(s0 * s1) * I01;
                double IG = 0.5 * Math.Log(m00 * m11 - m01 * m01);
                // kappa_eff: largest eigenvalue of (T I)^{1/2} Cov (T I)^{1/2}: use eigenvalues of T*I*Cov (same spectrum)
                double a00 = mT * (I00 * c00 + I01 * c01), a01 = mT * (I00 * c01 + I01 * c11), a10 = mT * (I01 * c00 + I11 * c01), a11 = mT * (I01 * c01 + I11 * c11);
                double tr = a00 + a11, dt = a00 * a11 - a01 * a10, kap = tr / 2 + Math.Sqrt(Math.Max(0, tr * tr / 4 - dt));
                double ldI = Math.Log(I00 * I11 - I01 * I01);
                Console.WriteLine(String.Format(ci, "{0,9:F1} {1,8:F3} {2,6:F3} {3,8:F3} {4,8:F3} {5,8:F3} {6,8:F3} {7,10:E3} {8,10:E3} {9,10:E3} {10,7:F3} {11,7:F3} {12,8:F3} {13,9:F0}",
                    mT, mR, se, L, mR - L, IG, mR - IG, c00, c11, c01, kap, ldI, Math.Log(detC * mT * mT), all.Average(x => x.n[g])));
            }
        }
        else
        {   // two-point: constant predictor theta' = theta + delta u, stationary start; directions: ridge-like u=(1,-w_r), top-like u=(0,1), mixed
            double al_s = (1 - eps) * beta; int warm = (int)(20.0 / (eps * eps)); double delta = double.Parse(get("delta", "0.1"), ci);
            // I from long run
            Console.WriteLine("# twopoint eps=" + eps + " mu=" + mu + " beta=" + beta + " delta=" + delta + "  T=" + tmax);
            Console.WriteLine("# direction  E_R/T  (u'Iu delta^2/2)  ratio  c-/c+ range of lam'/lam  1/c~+  1/c~-^2");
            var dirs = new List<double[]> { new double[] { 1, 0 }, new double[] { 0, 1 }, new double[] { 1, -eps * beta / mu }, new double[] { 1, 1 } };
            var sums = new double[dirs.Count][]; var I3 = new double[reps][]; var rmin = new double[dirs.Count]; var rmax = new double[dirs.Count]; for (int k = 0; k < dirs.Count; k++) { rmin[k] = double.MaxValue; rmax[k] = 0; }
            var Rk = new double[reps][];
            Parallel.For(0, reps, new ParallelOptions { MaxDegreeOfParallelism = 8 }, r =>
            {
                var rng = new Random(seed0 + r); double s = 0; for (int i = 0; i < warm; i++) { double ww = Gap(rng, mu, al_s, beta, s); s = s * Math.Exp(-beta * ww) + 1; }
                var Rr = new double[dirs.Count]; double Tt = 0, j0 = 0, j1 = 0, j2 = 0; var rmn = new double[dirs.Count]; var rmx = new double[dirs.Count]; for (int k = 0; k < dirs.Count; k++) { rmn[k] = double.MaxValue; }
                while (Tt < tmax)
                {
                    double wg = Gap(rng, mu, al_s, beta, s); double e = Math.Exp(-beta * wg), b1 = s * (1 - e) / beta;
                    for (int k = 0; k < dirs.Count; k++)
                    {
                        double mt = mu + delta * dirs[k][0], at = al_s + delta * dirs[k][1]; if (at >= beta) at = beta * 0.999;
                        double lg = 0; int ns = 4; double h = wg / ns;
                        for (int j = 0; j < ns; j++) for (int q = 0; q < 8; q++) { double u = h * (j + 0.5 * (gx[q] + 1)), su = s * Math.Exp(-beta * u); double ls = mu + al_s * su, lt = mt + at * su; lg += 0.5 * h * gw[q] * ls * Math.Log(lt / ls); double rt = lt / ls; if (rt < rmn[k]) rmn[k] = rt; if (rt > rmx[k]) rmx[k] = rt; }
                        Rr[k] += (mt - mu) * wg + (at - al_s) * b1 - lg;
                    }
                    double i0, i1, i2; FI(mu, al_s, beta, s, wg, out i0, out i1, out i2); j0 += i0; j1 += i1; j2 += i2; s = s * e + 1; Tt += wg;
                }
                Rk[r] = Rr; I3[r] = new double[] { j0 / Tt, j1 / Tt, j2 / Tt, Tt }; lock (rmin) { for (int k = 0; k < dirs.Count; k++) { rmin[k] = Math.Min(rmin[k], rmn[k]); rmax[k] = Math.Max(rmax[k], rmx[k]); } }
            });
            double I00 = I3.Average(x => x[0]), I01 = I3.Average(x => x[1]), I11 = I3.Average(x => x[2]), Tm = I3.Average(x => x[3]);
            Console.WriteLine("# Ibar = [" + I00.ToString("E4", ci) + ", " + I01.ToString("E4", ci) + "; ., " + I11.ToString("E4", ci) + "]  det=" + (I00 * I11 - I01 * I01).ToString("E4", ci));
            for (int k = 0; k < dirs.Count; k++)
            {
                double u0 = dirs[k][0] * delta, u1 = dirs[k][1] * delta; double qf = u0 * u0 * I00 + 2 * u0 * u1 * I01 + u1 * u1 * I11;
                double mR = Rk.Average(x => x[k]) / Tm, se = Math.Sqrt(Rk.Sum(x => (x[k] / Tm - mR) * (x[k] / Tm - mR)) / Math.Max(1, reps - 1) / reps);
                double cp = Math.Max(1, rmax[k]), cm = Math.Min(1, rmin[k]);
                Console.WriteLine(String.Format(ci, "u=({0:F2},{1:F3}) {2,10:E4}(+-{3:E1}) {4,10:E4} {5,7:F3}   [{6:F3},{7:F3}]  {8:F3} {9:F3}", dirs[k][0], dirs[k][1], mR, se, qf / 2, mR / (qf / 2), rmin[k], rmax[k], 1 / cp, 1 / (cm * cm)));
            }
        }
    }
}
