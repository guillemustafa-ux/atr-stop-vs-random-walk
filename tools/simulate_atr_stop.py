"""Sanity simulation of the ATR stops vs a random walk test on synthetic bars.

Standard library only. Run with:  python -I tools/simulate_atr_stop.py

Every series is built bar by bar from SUBSTEPS intrabar steps in log price: open = previous close, high and low
= the highest and the lowest point of the sub-path (the open included), close = its last point. Four series:

  (a) a driftless random walk in log price (geometric Brownian motion sampled at SUBSTEPS points per bar), the
      no-bias control: the measured hit rates should sit within a few points of the formula. The sub-path is
      watched at discrete points only, so it misses touches between them and the measured rates fall slightly
      below the formula; the Broadie-Glasserman-Kou continuity correction predicts by how much;
  (b) the same walk with FEW_SUBSTEPS points per bar, where that downward bias is large enough to see;
  (c) a mean-reverting series (Ornstein-Uhlenbeck in log price): the pull toward the mean caps how far price
      strays in N bars, so the wide stops are hit far less often than the formula says and the gap widens with
      k; the tightest stop, decided in the first few bars, stays near the formula;
  (d) a fat-tailed series (rare jumps on top of a smaller diffusion, same variance per bar): sigma is inflated by
      the jumps and the ATR much less, so the ATR / sigma ratio drops and the stops are hit less than predicted.

For every bar the simulation does what the Pine script does on a confirmed bar, in the same order: the open event
takes the bar's low and high into its lowest low and highest high; when its N bars are up it resolves (long stop
entry - k * ATR hit when the lowest low reaches it, short stop entry + k * ATR hit when the highest high reaches
it, and the random-walk prediction 2 * Phi(-d / (sigma * sqrt(N))) for each side, d the log distance from the
entry to the stop); then, with no event open and the ATR and the volatility window ready, a new event starts at
the close. The ATR is ta.atr: ta.rma of the true range, seeded with the simple average of the first ATR_LEN
values. Sigma is ta.stdev(log return, VOL_LEN, false): the sample standard deviation of the VOL_LEN 1-bar log
returns that end at the entry bar. Phi is the same Abramowitz and Stegun 7.1.26 approximation as in the script,
checked here against math.erf.

Printed per series, with the labels of the Pine table, then PASS or FAIL for every check, then the expected sign
of each result and what was observed.
"""

import math
import random
import statistics

SEED = 11
ATR_LEN = 14
VOL_LEN = 100  # V: bars of 1-bar log returns behind sigma
HORIZON = 20  # N: bars after the entry on which the stops are checked
KS = (1.0, 1.5, 2.0, 3.0)  # stop multiples of the ATR, ascending
BARS = 50000
SIGMA = 0.005  # standard deviation of a bar's log return
SUBSTEPS = 64  # intrabar points per bar
FEW_SUBSTEPS = 4
OU_THETA = 0.3  # mean reversion per bar of the Ornstein-Uhlenbeck series
JUMP_RATE = 0.02  # jumps per bar of the fat-tailed series
JUMP_SIZE = 6.0 * SIGMA  # log size of a jump, either sign
BGK_BETA = 0.5826  # -zeta(1/2) / sqrt(2 pi), the Broadie-Glasserman-Kou continuity correction
BROWNIAN_RANGE = math.sqrt(8.0 / math.pi)  # expected range of Brownian motion over one bar, in sigmas: 1.596
PHI_MAX_ERROR = 7.5e-8  # half the 1.5e-7 bound that Abramowitz and Stegun give for 7.1.26 on erf
TOLERANCE_PTS = 3.0  # "within a few points" of the formula
DISCRETE_TOLERANCE_PTS = 1.5  # about three standard errors of a pooled rate at 2,500 events

checks = []  # (passed, text)


def check(passed, text):
    checks.append((bool(passed), text))


def phi(x):
    """Same arithmetic as phi() in the Pine script: A&S 7.1.26 for erf, then Phi(x) = (1 + erf(x / sqrt 2)) / 2."""
    z = abs(x) / math.sqrt(2.0)
    t = 1.0 / (1.0 + 0.3275911 * z)
    erf_z = 1.0 - t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429)))) * math.exp(-z * z)
    return 0.5 * (1.0 + erf_z) if x >= 0 else 0.5 * (1.0 - erf_z)


def phi_exact(x):
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def hit_prob(d, s_n):
    """Chance that driftless Brownian motion touches a level at log distance d > 0 within N bars, s_n = sigma * sqrt(N)."""
    return 2.0 * phi(-d / s_n)


# ───────────────────────── Series ─────────────────────────
def bars_from(steps_of_bar):
    """Turns a generator of per-bar log sub-paths (open point first) into (open, high, low, close) bars."""
    for path in steps_of_bar:
        o, h, l, c = math.exp(path[0]), math.exp(max(path)), math.exp(min(path)), math.exp(path[-1])
        yield o, h, l, c


def walk(rng, substeps):
    x = math.log(100.0)
    s = SIGMA / math.sqrt(substeps)
    for _ in range(BARS):
        path = [x]
        for _ in range(substeps):
            x += s * rng.gauss(0.0, 1.0)
            path.append(x)
        yield path


def mean_reverting(rng, substeps):
    mu = x = math.log(100.0)
    s = SIGMA / math.sqrt(substeps)
    pull = OU_THETA / substeps
    for _ in range(BARS):
        path = [x]
        for _ in range(substeps):
            x += -pull * (x - mu) + s * rng.gauss(0.0, 1.0)
            path.append(x)
        yield path


def fat_tailed(rng, substeps):
    x = math.log(100.0)
    s = math.sqrt(SIGMA ** 2 - JUMP_RATE * JUMP_SIZE ** 2) / math.sqrt(substeps)  # same variance per bar as (a)
    for _ in range(BARS):
        jump_at = rng.randrange(substeps) if rng.random() < JUMP_RATE else -1
        sign = 1.0 if rng.random() < 0.5 else -1.0
        path = [x]
        for j in range(substeps):
            x += s * rng.gauss(0.0, 1.0) + (sign * JUMP_SIZE if j == jump_at else 0.0)
            path.append(x)
        yield path


# ───────────────────────── Counting, as in the Pine script ─────────────────────────
def measure(bars, substeps):
    """Runs the event ledger over the bars; returns per-k counts and the volatility summary."""
    k_n = len(KS)
    hit_l, hit_s, hit_both = [0] * k_n, [0] * k_n, [0] * k_n
    pred_l, pred_s, pred_disc = [0.0] * k_n, [0.0] * k_n, [0.0] * k_n
    n_ev = 0
    sum_atr_pct = sum_sig_pct = 0.0
    tr_seed = []  # first ATR_LEN true ranges, for the SMA seed of ta.rma
    atr = None
    rets = []  # 1-bar log returns, None on the first bar
    prev_close = None
    ev = None  # open event: [bar, entry, atr, sigma, lowest low, highest high]
    i = -1
    for i, (o, h, l, c) in enumerate(bars):
        tr = h - l if prev_close is None else max(h - l, abs(h - prev_close), abs(l - prev_close))
        if atr is None:
            tr_seed.append(tr)
            if len(tr_seed) == ATR_LEN:
                atr = sum(tr_seed) / ATR_LEN
        else:
            atr = (tr + (ATR_LEN - 1) * atr) / ATR_LEN
        rets.append(None if prev_close is None else math.log(c / prev_close))
        prev_close = c
        if ev is not None:
            ev[4] = min(ev[4], l)
            ev[5] = max(ev[5], h)
            if i - ev[0] >= HORIZON:
                _, e, a, sg, lo, hi = ev
                s_n = sg * math.sqrt(HORIZON)
                shift = BGK_BETA * sg / math.sqrt(substeps)  # watched at `substeps` points per bar
                for j, k in enumerate(KS):
                    stop_l, stop_s = e - k * a, e + k * a
                    hl, hs = lo <= stop_l, hi >= stop_s
                    hit_l[j] += hl
                    hit_s[j] += hs
                    hit_both[j] += hl and hs
                    d_l = math.log(e / stop_l) if stop_l > 0 else math.inf
                    d_s = math.log(stop_s / e)
                    pred_l[j] += hit_prob(d_l, s_n) if stop_l > 0 else 0.0
                    pred_s[j] += hit_prob(d_s, s_n)
                    pred_disc[j] += (hit_prob(d_l + shift, s_n) if stop_l > 0 else 0.0) + hit_prob(d_s + shift, s_n)
                n_ev += 1
                sum_atr_pct += a / e * 100
                sum_sig_pct += sg * 100
                ev = None
        if ev is None and i >= VOL_LEN and atr is not None:
            sg = statistics.stdev(rets[i - VOL_LEN + 1:i + 1])
            if atr > 0 and sg > 0:
                ev = [i, c, atr, sg, c, c]
    rows = []
    for j, k in enumerate(KS):
        n = n_ev
        p_l, p_s = hit_l[j] / n, hit_s[j] / n
        p = (hit_l[j] + hit_s[j]) / (2 * n)
        # The event is the unit: y = (long hit + short hit) / 2 per event, so the interval keeps the dependence
        # between the two sides of one event. Sum of y^2 = (long hits + short hits + 2 * both) / 4.
        var = max((hit_l[j] + hit_s[j] + 2 * hit_both[j]) / (4 * n) - p * p, 0.0)
        half = 1.96 * math.sqrt(var / n)
        rw = (pred_l[j] + pred_s[j]) / (2 * n)
        rw_disc = pred_disc[j] / (2 * n)
        rows.append(dict(k=k, n=n, hl=hit_l[j], hs=hit_s[j], p_l=p_l, p_s=p_s, p=p, half=half, rw=rw, rw_disc=rw_disc))
    atr_pct, sig_pct = sum_atr_pct / n_ev, sum_sig_pct / n_ev
    return rows, dict(n_ev=n_ev, bars=i + 1, atr_pct=atr_pct, sig_pct=sig_pct, ratio=atr_pct / sig_pct)


def verdict(r):
    """Same rule as the Pine script: the random-walk rate against the pooled rate's 95% interval."""
    if r["n"] == 0:
        return "not enough events yet"
    v = ("hit more than a random walk" if r["rw"] < r["p"] - r["half"] else
         "hit less than a random walk" if r["rw"] > r["p"] + r["half"] else "like a random walk")
    return v + (" (small sample)" if r["n"] < 30 else "")


def pct(p):
    """Same template as fmtRate() in the Pine script ("0.0" and a percent sign)."""
    return f"{100.0 * p:.1f}%"


def pts(v):
    """Same template as fmtPts() in the Pine script: signed points, "0.0 pts" when it rounds to zero."""
    return "0.0 pts" if abs(v) < 0.05 else f"{v:+.1f} pts"


def report(name, rows, summ, substeps):
    print(f"{name}  ({summ['bars']} bars, {substeps} sub-steps per bar, seed {SEED}, ATR {ATR_LEN}, sigma over "
          f"{VOL_LEN} bars, horizon {HORIZON})")
    print(f"  {'Stop':<10}{'events':>7}{'long hit':>10}{'short hit':>11}{'pooled hit':>12}   {'95% CI':<16}"
          f"{'random walk':>11}{'difference':>12}   verdict   [discrete-watch formula]")
    for r in rows:
        ci = f"{pct(max(r['p'] - r['half'], 0.0))} - {pct(min(r['p'] + r['half'], 1.0))}"
        label = f"{r['k']:g} x ATR"
        print(f"  {label:<10}{r['n']:>7}{pct(r['p_l']):>10}{pct(r['p_s']):>11}{pct(r['p']):>12}   {ci:<16}"
              f"{pct(r['rw']):>11}{pts(100 * (r['p'] - r['rw'])):>12}   {verdict(r)}   [{pct(r['rw_disc'])}]")
    print(f"  Volatility  ATR {summ['atr_pct']:.2f}% of price   sigma per bar {summ['sig_pct']:.2f}%   "
          f"ATR / sigma {summ['ratio']:.2f}   (Brownian bar range {BROWNIAN_RANGE:.2f})   {summ['bars']} confirmed bars")
    print()


def identities(name, rows, summ):
    """The relations the Pine table must satisfy, checked on the simulated counts."""
    ok = all(r["n"] == summ["n_ev"] for r in rows)
    ok &= all(0 <= r["hl"] <= r["n"] and 0 <= r["hs"] <= r["n"] for r in rows)
    ok &= all(abs(r["p"] - (r["p_l"] + r["p_s"]) / 2) < 1e-12 for r in rows)
    ok &= all(a["hl"] >= b["hl"] and a["hs"] >= b["hs"] and a["rw"] >= b["rw"] for a, b in zip(rows, rows[1:]))
    ok &= abs(summ["ratio"] - summ["atr_pct"] / summ["sig_pct"]) < 1e-12
    check(ok, f"{name}: identities (same events in every row; hits <= events; pooled = (long + short) / 2; "
              "hits and the random-walk rate never rise with k; ATR / sigma = ATR % / sigma %)")


def main():
    # Phi against math.erf on a grid of 32,001 points over [-8, 8].
    grid = [i / 2000.0 for i in range(-16000, 16001)]
    worst = max(abs(phi(x) - phi_exact(x)) for x in grid)
    check(worst <= PHI_MAX_ERROR, f"Phi approximation: max |Phi_A&S - Phi_erf| = {worst:.2e} over [-8, 8] "
                                  f"(bound {PHI_MAX_ERROR:.1e})")

    rng = random.Random(SEED)
    gbm_rows, gbm = measure(bars_from(walk(rng, SUBSTEPS)), SUBSTEPS)
    report("Random walk", gbm_rows, gbm, SUBSTEPS)
    rng = random.Random(SEED)
    few_rows, few = measure(bars_from(walk(rng, FEW_SUBSTEPS)), FEW_SUBSTEPS)
    report("Random walk, few sub-steps", few_rows, few, FEW_SUBSTEPS)
    rng = random.Random(SEED)
    ou_rows, ou = measure(bars_from(mean_reverting(rng, SUBSTEPS)), SUBSTEPS)
    report("Mean-reverting", ou_rows, ou, SUBSTEPS)
    rng = random.Random(SEED)
    ft_rows, ft = measure(bars_from(fat_tailed(rng, SUBSTEPS)), SUBSTEPS)
    report("Fat-tailed", ft_rows, ft, SUBSTEPS)

    for name, rows, summ in (("Random walk", gbm_rows, gbm), ("Few sub-steps", few_rows, few),
                             ("Mean-reverting", ou_rows, ou), ("Fat-tailed", ft_rows, ft)):
        identities(name, rows, summ)

    gbm_diffs = [100 * (r["p"] - r["rw"]) for r in gbm_rows]
    check(all(abs(d) <= TOLERANCE_PTS for d in gbm_diffs),
          f"Random walk: measured within {TOLERANCE_PTS:g} points of the formula on every k "
          f"(differences {', '.join(f'{d:+.1f}' for d in gbm_diffs)})")
    gbm_disc = [100 * (r["p"] - r["rw_disc"]) for r in gbm_rows]
    check(all(abs(d) <= DISCRETE_TOLERANCE_PTS for d in gbm_disc),
          f"Random walk: measured within {DISCRETE_TOLERANCE_PTS:g} points of the discrete-watch formula on every k "
          f"({', '.join(f'{d:+.1f}' for d in gbm_disc)})")
    few_diffs = [100 * (r["p"] - r["rw"]) for r in few_rows]
    few_disc = [100 * (r["p"] - r["rw_disc"]) for r in few_rows]
    check(all(d < -2.0 for d in few_diffs) and all(abs(d) <= DISCRETE_TOLERANCE_PTS for d in few_disc),
          f"Few sub-steps: measured more than 2 points below the formula on every k "
          f"({', '.join(f'{d:+.1f}' for d in few_diffs)}), and within {DISCRETE_TOLERANCE_PTS:g} of the discrete-watch "
          f"formula ({', '.join(f'{d:+.1f}' for d in few_disc)})")
    expected_ratio = BROWNIAN_RANGE - 2 * BGK_BETA / math.sqrt(SUBSTEPS)
    check(abs(gbm["ratio"] - expected_ratio) <= 0.03,
          f"Random walk: ATR / sigma {gbm['ratio']:.3f} within 0.03 of {expected_ratio:.3f}, the Brownian range "
          f"{BROWNIAN_RANGE:.3f} less the discrete-watch shortfall at {SUBSTEPS} sub-steps")
    ou_diffs = [100 * (r["p"] - r["rw"]) for r in ou_rows]
    check(all(a > b for a, b in zip(ou_diffs, ou_diffs[1:])) and ou_diffs[-1] < -10.0
          and verdict(ou_rows[-1]).startswith("hit less"),
          f"Mean-reverting: the difference falls as k grows ({', '.join(f'{d:+.1f}' for d in ou_diffs)}) and the "
          f"widest stop reads 'hit less than a random walk', more than 10 points below")
    ft_diffs = [100 * (r["p"] - r["rw"]) for r in ft_rows]
    check(all(verdict(r).startswith("hit less") for r in ft_rows) and ft["ratio"] < gbm["ratio"] - 0.2,
          f"Fat-tailed: 'hit less than a random walk' on every k ({', '.join(f'{d:+.1f}' for d in ft_diffs)}) "
          f"and ATR / sigma {ft['ratio']:.2f} well below the random walk's {gbm['ratio']:.2f}")

    for passed, text in checks:
        print(f"{'PASS' if passed else 'FAIL'}  {text}")
    print()
    print("Expected sign: on the random walk the measured rates sit a little below the formula, by about what the "
          "discrete-watch correction predicts. "
          f"Observed: mean difference {sum(gbm_diffs) / len(gbm_diffs):+.1f} pts at {SUBSTEPS} sub-steps, "
          f"{sum(few_diffs) / len(few_diffs):+.1f} pts at {FEW_SUBSTEPS}.")
    print("Expected sign: the mean-reverting series hits the wide stops less than the formula says, more so as k "
          "grows; the fat-tailed series hits every stop less, with an ATR / sigma ratio well below the random walk's. "
          f"Observed: {ou_diffs[0]:+.1f} pts at k = {KS[0]:g} down to {ou_diffs[-1]:+.1f} pts at k = {KS[-1]:g}; "
          f"fat-tailed mean difference {sum(ft_diffs) / len(ft_diffs):+.1f} pts, ratio {ft['ratio']:.2f}.")
    failed = sum(not p for p, _ in checks)
    print(f"\n{len(checks) - failed} of {len(checks)} checks PASS.")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
