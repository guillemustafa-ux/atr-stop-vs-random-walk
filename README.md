# ATR Stops vs a Random Walk

Open-source TradingView indicator (Pine Script v6) that places a long and a short stop at four ATR multiples from the close at regular, non-overlapping entries, counts how often the next N bars reach each one, and compares every hit rate with what a driftless random walk with the same volatility would give.

## What it measures

A stop k ATRs away is meant to sit outside the noise. How often noise alone reaches it has a closed form: for Brownian motion in log price with volatility sigma per bar, the chance of touching a level at log distance `d` within N bars is `2 * Phi(-d / (sigma * sqrt(N)))`, by the reflection principle. The indicator measures the real hit rate on the chart and puts that prediction next to it:

1. **Starts events.** Once the ATR and the volatility window are ready, an event starts at the close of a confirmed bar and lasts N bars (default 20). The next event starts on the bar where it resolves (or, if the ATR or sigma is zero there, on the next bar where both are above zero), so the entries are at least N bars apart and the windows never overlap. Each event keeps its entry (the close), the ATR of the entry bar, `ta.atr(atrLen)` (default 14), and sigma, `ta.stdev(math.log(close / close[1]), volLen, false)`: the sample standard deviation of the `volLen` (default 100) one-bar log returns that end at the entry bar. Both are known at the close of the entry bar.
2. **Checks the stops.** Four multiples `k` (defaults 1, 1.5, 2 and 3). The long stop is `entry - k * ATR` and is hit when a low of the N bars after the entry is at or below it; the short stop is `entry + k * ATR` and is hit when a high is at or above it. Both sides are measured at every event, so each event gives one long and one short outcome per `k`.
3. **Predicts.** Per event, side and `k`: `d = ln(entry / long stop)` or `d = ln(short stop / entry)`, and the prediction `2 * Phi(-d / (sigma * sqrt(N)))` with the event's own sigma. A long stop at or below zero cannot be reached by a positive price and gets 0. The random-walk rate of a row is the mean of the long and the short predictions over all events.
4. **Reports.** On the last bar, per `k`: events, long, short and pooled hit rates, the pooled rate's 95% interval, the random-walk rate, the difference in points and a verdict; then the volatility row and the last event.

Pine has no normal distribution function, so `Phi` is built from erf with formula 7.1.26 of Abramowitz and Stegun, *Handbook of Mathematical Functions* (1964): for `x >= 0`, `erf(x) = 1 - (a1 t + a2 t^2 + a3 t^3 + a4 t^4 + a5 t^5) exp(-x^2)`, `t = 1 / (1 + p x)`, `p = 0.3275911`, `a1..a5 = 0.254829592, -0.284496736, 1.421413741, -1.453152027, 1.061405429`, with an absolute error of at most 1.5e-7. With `Phi(x) = (1 + erf(x / sqrt(2))) / 2` and `erf(-x) = -erf(x)`, the error on `Phi` is at most 7.5e-8; `tools/simulate_atr_stop.py` runs the same arithmetic against `math.erf` on 32,001 points over [-8, 8] and measured a maximum of 6.97e-8.

Why a volatility row: the prediction takes sigma from close-to-close returns, the stops take the ATR from bar ranges. For a Brownian bar the expected range is `sqrt(8 / pi) * sigma`, about 1.60 sigma: by the reflection principle the maximum over the bar is distributed as `|Z| * sigma`, whose mean is `sqrt(2 / pi) * sigma`, and the minimum mirrors it. So on a random-walk chart the ATR is about 1.6 sigma and a stop k ATRs away is about 1.6 k sigmas away. The row shows the mean ATR as a percentage of the entry, the mean sigma and the ratio of the two: a ratio far from 1.6 says the bars do not move like a random walk inside the bar, whatever the stops do, and the verdicts reflect that too. In the simulation a fat-tailed series (rare jumps) reads 0.92 and is hit 12 to 16 points less than the formula at every `k`.

## Counting rules

The counters update only on confirmed bars (`barstate.isconfirmed`), in this order:

1. **Window.** The open event takes the bar's low and high into its lowest low and highest high. Both start at the entry close, which is above every long stop and below every short stop, so the start value never decides a hit.
2. **Resolution.** When N bars have passed since the entry bar, every stop is compared with the lowest low or the highest high, the predictions are added, and the event closes.
3. **New event.** With no event open, `bar_index >= volLen`, and the ATR and sigma available and above zero, a new event starts at this close.

Consequences:

- The realtime bar is never in the counts until it confirms. The event still open on the last bar is not counted; the "Last event" row names it.
- The entry bar's own high and low are not in its window: the window is the N bars after it, and the next entry bar is the last bar of the window before it.
- Entries are spaced from the first bar where the ATR and the window are ready, so the schedule depends on the first loaded bar; a different history depth shifts every entry.
- A stop is hit on a touch. Which of the long and the short stop came first is not tracked: one event can hit both.
- No `request.*` call, no look-ahead, no recalculation of a count after its bar has closed.

## Statistics

- For a row with `n` events, `hL` long hits, `hS` short hits and `hB` events where both were hit: long rate `hL / n`, short rate `hS / n`, pooled rate `p = (hL + hS) / (2n)`.
- 95% interval, with the event as the unit: per event `y = (long hit + short hit) / 2` is 0, 0.5 or 1 and its mean is `p`; `Var(y) = (hL + hS + 2 hB) / (4n) - p^2`; interval `p ± 1.96 * sqrt(Var(y) / n)`, clipped to [0, 1] for display. The two sides of one event share the same bars and are not independent; when they happen to be independent this is the binomial interval on `2n` trials.
- Random-walk rate: the mean of the `2n` predictions. Difference: `(p - random-walk rate) * 100`, in points.
- Verdict: "like a random walk" when the random-walk rate lies inside the interval; "hit more than a random walk" when it lies below the interval's lower bound; "hit less than a random walk" when it lies above the upper bound; "not enough events yet" while no event has resolved. Under 30 events the verdict adds "(small sample)".
- The random-walk rate is treated as a fixed number. The sigma behind it is an estimate with a relative error of about `1 / sqrt(2 V)`, 7% at V = 100, which the interval does not include.
- Volatility row: mean of `ATR / entry * 100` over the events, mean of `sigma * 100`, and the ratio of the two means.

## Settings

| Input | Default | Meaning |
| --- | --- | --- |
| ATR length | 14 | Length of `ta.atr(length)`, read on the entry bar. |
| Stop 1 (x ATR) | 1.0 | Multiple `k` of the first row. Keep the four stops in ascending order so the rows read from the closest stop to the widest. |
| Stop 2 (x ATR) | 1.5 | Multiple `k` of the second row. |
| Stop 3 (x ATR) | 2.0 | Multiple `k` of the third row. |
| Stop 4 (x ATR) | 3.0 | Multiple `k` of the fourth row. |
| Horizon N (bars) | 20 | Bars after the entry on which the stops are checked, 1 to 500; also the spacing of the entries. |
| Volatility window V (bars) | 100 | One-bar log returns behind sigma, 10 to 500. |
| Mark entries | on | A small circle below each entry bar. |
| Table position | Bottom right | Four corners. |
| Table text size | Small | Tiny, Small or Normal. |

## How to read the table

| Row | Meaning |
| --- | --- |
| Header | Symbol, timeframe, ATR length, volatility window, horizon. |
| Column titles | Stop, events, long hit, short hit, pooled hit, 95% CI, random walk, difference, verdict. |
| One row per stop | The four `k` rows in input order. |
| Volatility | Mean ATR as % of price, mean sigma per bar, their ratio, the Brownian reference 1.60, confirmed bars and events. |
| Last event | Entry time, close, ATR, the lowest low and the highest high of its N bars in ATRs from the entry ("low -1.84 ATR" means the long stops up to 1.84 x ATR were hit), and the open event that is not counted yet. |

The footer reminds that the numbers are counts on confirmed bars and names the timezone of the dates (`syminfo.timezone`). Rates are shown with one decimal.

## Files

- `atr_stop_vs_random_walk.pine`: the indicator.
- `TRADINGVIEW_PUBLICATION.md`: title, description text and checklist for publishing.
- `tools/simulate_atr_stop.py`: standard-library Python simulation (`python -I tools/simulate_atr_stop.py`) of the same ledger (ATR as `ta.atr`, sigma as `ta.stdev(..., false)`, the same `Phi` arithmetic, the same order on a bar) on four seeded series of 50,000 bars built from intrabar sub-steps in log price, where the open is the previous close and the high and the low are the extremes of the sub-path: a driftless random walk at 64 sub-steps per bar, the same walk at 4 sub-steps per bar, a mean-reverting series and a fat-tailed series. It prints the table rows of each and PASS or FAIL for eleven checks: the `Phi` approximation against `math.erf`, the table identities on every series, the random walk within 3 points of the formula and within 1.5 points of the discrete-watch formula, the larger shortfall at 4 sub-steps, the random walk's ATR / sigma ratio, and the departure of the mean-reverting and the fat-tailed series.
- `LICENSE`: MIT.

## Check it yourself

1. Add the indicator to BINANCE:BTCUSDT on 1 hour from the Pine Editor and confirm it compiles without errors or warnings.
2. Check the identities: no rate above 100%; the same number of events in every stop row; the pooled rate is the mean of the long and the short rates up to rounding; with the default ascending stops the long, the short and the random-walk rates never rise from one row to the next (a stop further away is reached only if the closer one was); "difference" equals pooled minus random walk up to rounding; ATR / sigma equals the ATR percentage over the sigma percentage up to rounding.
3. Set Stop 1 to 0.1: the first row's pooled and random-walk rates must read near 100% (about 97% at N = 20 with a ratio near 1.6). Near 100% the normal interval is unreliable (zero width when every event hits), so the verdict on that row means nothing. Set it back to 1.
4. Set the horizon to 1: every row must be small (the random walk gives about 11% at 1 ATR with a ratio of 1.6, and less further out), and there are about 20 times the events. Set it back to 20.
5. Hand check the last event: find the circle at its entry time; over the next 20 bars, the lowest low minus the entry close, divided by the ATR in the row, must equal the "low" cell (0.00 when no low went below the entry), and the highest high minus the entry close, divided by the ATR, the "high" cell (0.00 when no high went above it).
6. Compare with `tools/simulate_atr_stop.py`. The symbol differs, so the numbers will too; what carries over is how to read them. Output on 2026-10-14:

```text
Random walk  (50000 bars, 64 sub-steps per bar, seed 11, ATR 14, sigma over 100 bars, horizon 20)
  Stop       events  long hit  short hit  pooled hit   95% CI          random walk  difference   verdict   [discrete-watch formula]
  1 x ATR      2494     72.9%      73.7%       73.3%   72.3% - 74.3%         74.3%    -1.0 pts   like a random walk   [73.0%]
  1.5 x ATR    2494     61.3%      61.7%       61.5%   60.6% - 62.4%         62.3%    -0.7 pts   like a random walk   [61.1%]
  2 x ATR      2494     49.4%      50.7%       50.1%   49.3% - 50.9%         51.2%    -1.1 pts   hit less than a random walk   [50.1%]
  3 x ATR      2494     31.1%      33.5%       32.3%   31.3% - 33.2%         32.6%    -0.3 pts   like a random walk   [31.8%]
  Volatility  ATR 0.73% of price   sigma per bar 0.50%   ATR / sigma 1.46   (Brownian bar range 1.60)   50000 confirmed bars

Random walk, few sub-steps  (50000 bars, 4 sub-steps per bar, seed 11, ATR 14, sigma over 100 bars, horizon 20)
  Stop       events  long hit  short hit  pooled hit   95% CI          random walk  difference   verdict   [discrete-watch formula]
  1 x ATR      2494     74.1%      76.2%       75.1%   74.1% - 76.1%         80.2%    -5.1 pts   hit less than a random walk   [75.2%]
  1.5 x ATR    2494     64.3%      67.6%       66.0%   65.0% - 66.9%         70.7%    -4.7 pts   hit less than a random walk   [65.9%]
  2 x ATR      2494     55.5%      59.2%       57.3%   56.5% - 58.2%         61.6%    -4.3 pts   hit less than a random walk   [57.1%]
  3 x ATR      2494     41.2%      43.5%       42.3%   41.5% - 43.2%         45.3%    -2.9 pts   hit less than a random walk   [41.5%]
  Volatility  ATR 0.56% of price   sigma per bar 0.50%   ATR / sigma 1.12   (Brownian bar range 1.60)   50000 confirmed bars

Mean-reverting  (50000 bars, 64 sub-steps per bar, seed 11, ATR 14, sigma over 100 bars, horizon 20)
  Stop       events  long hit  short hit  pooled hit   95% CI          random walk  difference   verdict   [discrete-watch formula]
  1 x ATR      2494     76.2%      77.3%       76.7%   75.8% - 77.7%         73.2%    +3.5 pts   hit more than a random walk   [72.0%]
  1.5 x ATR    2494     57.0%      57.7%       57.4%   56.4% - 58.3%         60.8%    -3.4 pts   hit less than a random walk   [59.6%]
  2 x ATR      2494     37.7%      38.2%       37.9%   36.9% - 38.9%         49.4%   -11.5 pts   hit less than a random walk   [48.4%]
  3 x ATR      2494      9.9%      11.0%       10.5%   9.7% - 11.3%          30.5%   -20.1 pts   hit less than a random walk   [29.8%]
  Volatility  ATR 0.71% of price   sigma per bar 0.47%   ATR / sigma 1.53   (Brownian bar range 1.60)   50000 confirmed bars

Fat-tailed  (50000 bars, 64 sub-steps per bar, seed 11, ATR 14, sigma over 100 bars, horizon 20)
  Stop       events  long hit  short hit  pooled hit   95% CI          random walk  difference   verdict   [discrete-watch formula]
  1 x ATR      2494     72.2%      69.6%       70.9%   69.9% - 71.9%         82.7%   -11.8 pts   hit less than a random walk   [81.4%]
  1.5 x ATR    2494     61.3%      58.7%       60.0%   59.1% - 60.9%         74.4%   -14.3 pts   hit less than a random walk   [73.1%]
  2 x ATR      2494     51.8%      49.6%       50.7%   49.8% - 51.7%         66.4%   -15.6 pts   hit less than a random walk   [65.2%]
  3 x ATR      2494     37.4%      35.4%       36.4%   35.4% - 37.5%         51.8%   -15.4 pts   hit less than a random walk   [50.8%]
  Volatility  ATR 0.44% of price   sigma per bar 0.48%   ATR / sigma 0.92   (Brownian bar range 1.60)   50000 confirmed bars

PASS  Phi approximation: max |Phi_A&S - Phi_erf| = 6.97e-08 over [-8, 8] (bound 7.5e-08)
PASS  Random walk: identities (same events in every row; hits <= events; pooled = (long + short) / 2; hits and the random-walk rate never rise with k; ATR / sigma = ATR % / sigma %)
PASS  Few sub-steps: identities (same events in every row; hits <= events; pooled = (long + short) / 2; hits and the random-walk rate never rise with k; ATR / sigma = ATR % / sigma %)
PASS  Mean-reverting: identities (same events in every row; hits <= events; pooled = (long + short) / 2; hits and the random-walk rate never rise with k; ATR / sigma = ATR % / sigma %)
PASS  Fat-tailed: identities (same events in every row; hits <= events; pooled = (long + short) / 2; hits and the random-walk rate never rise with k; ATR / sigma = ATR % / sigma %)
PASS  Random walk: measured within 3 points of the formula on every k (differences -1.0, -0.7, -1.1, -0.3)
PASS  Random walk: measured within 1.5 points of the discrete-watch formula on every k (+0.3, +0.4, -0.1, +0.5)
PASS  Few sub-steps: measured more than 2 points below the formula on every k (-5.1, -4.7, -4.3, -2.9), and within 1.5 of the discrete-watch formula (-0.1, +0.0, +0.2, +0.9)
PASS  Random walk: ATR / sigma 1.464 within 0.03 of 1.450, the Brownian range 1.596 less the discrete-watch shortfall at 64 sub-steps
PASS  Mean-reverting: the difference falls as k grows (+3.5, -3.4, -11.5, -20.1) and the widest stop reads 'hit less than a random walk', more than 10 points below
PASS  Fat-tailed: 'hit less than a random walk' on every k (-11.8, -14.3, -15.6, -15.4) and ATR / sigma 0.92 well below the random walk's 1.46

Expected sign: on the random walk the measured rates sit a little below the formula, by about what the discrete-watch correction predicts. Observed: mean difference -0.8 pts at 64 sub-steps, -4.3 pts at 4.
Expected sign: the mean-reverting series hits the wide stops less than the formula says, more so as k grows; the fat-tailed series hits every stop less, with an ATR / sigma ratio well below the random walk's. Observed: +3.5 pts at k = 1 down to -20.1 pts at k = 3; fat-tailed mean difference -14.3 pts, ratio 0.92.

11 of 11 checks PASS.
```

The bracketed column exists only in the simulation. The synthetic high and low are the extremes of 64 (or 4) points per bar, not of a continuous path, so they miss touches between the points and the measured rates fall below the formula. Broadie, Glasserman and Kou (1997, "A continuity correction for discrete barrier options", *Mathematical Finance* 7(4)) give the size of that shortfall: a path watched at steps of `dt` bars hits a level at distance `d` about as often as a continuous one hits a level at `d + 0.5826 * sigma * sqrt(dt)` (0.5826 = `-zeta(1/2) / sqrt(2 pi)`). The bracketed rate is the formula with that shift; it lands within 1.5 points of the measured rate on every row of both random walks, and the same shortfall explains why the simulated ATR / sigma ratio is 1.46 at 64 sub-steps and not 1.60: the expected range of a bar watched at M points is about `(1.596 - 1.165 / sqrt(M)) * sigma`, about 1.45 at M = 64 (exact 1.456 by Spitzer's sum, `E[max] = sum over k = 1..M of E[S_k+] / k`).

Three things to notice in the output. At 64 sub-steps the random walk sits 0.3 to 1.1 points below the formula, and that is enough for one row to read "hit less than a random walk" on a sample of 2,494 events: the discrete-watch formula, inside the interval on that row, shows the shortfall is the sub-steps', not the stop's. On a real chart the high and the low are the extremes of the traded prices, so this shortfall is negligible when a bar holds many trades; on thin symbols it returns, and the bid-ask bounce pushes the other way. The mean-reverting series is hit more than the formula at 1 ATR and far less at 3 ATRs. Part of it can be checked: mean reversion lowers the close-to-close sigma below the intrabar diffusion's, `Var(1-bar return) = sigma^2 * (1 - e^-theta) / theta`, 0.864 sigma^2 at theta = 0.3, a sigma of 0.929 times the diffusion's (the simulation measured 0.932). The formula takes sigma from closes, so it under-predicts the tight stop: with the diffusion's sigma the 1 ATR prediction is 75.0% instead of 73.2%, 1.8 of the 3.5 points, and the rest is not explained by sigma alone. At the wide stops the pull caps how far price strays: the stationary spread is `sigma / sqrt(2 theta)` = 1.29 diffusion sigmas, the 3 ATR stop sits about 4.3 diffusion sigmas from the entry, and even with the diffusion's sigma the formula gives 34.1% against the 10.5% measured. The fat-tailed series is hit 12 to 16 points less at every stop, and its ATR / sigma ratio of 0.92 flags it before any stop is read.

Result on TradingView: not compiled yet. The script has to be compiled in the Pine Editor and read on BINANCE:BTCUSDT 1 hour before the numbers can be reported here.

## Limitations

- The null is a driftless random walk with a constant volatility per event. A trend, a change of volatility, mean reversion or fat tails can move the hit rates away from it; the verdict says that the chart differs from the walk, not why. Pooling both sides cancels most of a steady drift, but a strong trend in the loaded bars still moves the pooled rate.
- Sigma is estimated from V returns and has its own error, about 7% at V = 100, which the interval does not include.
- The ATR (an RMA of 14 bars by default) and sigma (100 bars by default) react to a change of volatility at different speeds, so around a volatility shift the stop distance and the prediction are measured on different windows.
- A hit is a touch of the high or the low. No fills, no slippage, no costs, and no order between the long and the short stop of one event.
- Highs and lows carry the bid-ask bounce one way (the top trades at the ask, the bottom at the bid), closes with either sign; on thin symbols or low timeframes it inflates the ATR against sigma and the hits at tight stops.
- The four rows share the same events, so they are not independent tests.
- Standard candles only: Heikin Ashi, Renko and other synthetic bars change the highs, the lows and the closes the counts read.
- Only bars loaded on the chart are counted, and the entries are spaced from the first ready bar, so history depth, symbol and timeframe change every count.
- It reports counts, not a trading rule. Not financial advice.
- The script has not been compiled outside TradingView; compile it in the Pine Editor before relying on it.

## Pine API notes

What the Pine Script v6 reference and user manual say about the built-ins the script relies on, read on 2026-10-14 through the context7 documentation index of the reference (`/websites/tradingview_pine-script-reference_v6`) and the user manual (`/websites/tradingview_pine-script-docs`):

- `ta.atr(length) → series float`, `length` a simple int: "returns the RMA of true range. True range is max(high - low, abs(high - close[1]), abs(low - close[1]))". The reference's equivalent code uses `high - low` when `high[1]` is na (the first bar) and then `ta.rma(trueRange, length)`. The script passes an input int.
- `ta.rma(source, length)`: the reference's equivalent code is `sum := na(sum[1]) ? ta.sma(src, length) : alpha * src + (1 - alpha) * nz(sum[1])` with `alpha = 1/length`, so the first value is the simple average of the first `length` values. The simulation seeds its ATR the same way.
- `ta.stdev(source, length, biased) → series float`: "If biased is true, function will calculate using a biased estimate of the entire population, if false - unbiased estimate of a sample"; the default is true and the script passes false, the sample estimate the simulation uses (`statistics.stdev`). Remarks: "na values in the source series are ignored; the function calculates on the length quantity of non-na values." The log return is na on the first bar, and an event needs `bar_index >= volLen`, so its window holds `volLen` real returns.
- `math.log(number)`: "Natural logarithm of any number > 0"; `math.exp(number)`: e raised to the power of `number`; `math.sqrt(number)`: the square root of `number >= 0`. The script calls `math.log` on ratios of positive prices only.
- `barstate.isconfirmed`: "Returns true if the script is calculating the last (closing) update of the current bar." All counters update inside it.
- `table.new(position, columns, rows, ...)`: "it is wise to use table.new() in conjunction with either the var keyword (so the table object is only created on the first bar) or in an if barstate.islast block"; the script uses `var`. `table.merge_cells(...)`: "An attempt to merge a cell that has already been merged will result in an error", which is why the merges run once behind the `merged` flag.
- `plotshape(series, title, style, location, color, offset, text, textcolor, editable, size, show_last, display, format, precision, force_overlay) → void`: styles include `shape.circle`, locations `location.belowbar`, sizes `size.tiny`; `display` takes `display.pane` among its values, which keeps the marker's value out of the status line; "Use plotshape() function in conjunction with 'overlay=true' indicator() parameter!"
- `str.tostring(value, format)`: "The formatting of float values will also round those values when necessary ... To display trailing zeros, use '0' instead of '#'", so `"0.###"` prints the stop multiples with up to three decimals and no trailing zeros; with `format.mintick` "the value will be rounded to the nearest number that can be divided by syminfo.mintick". `str.format_time(time, format, timezone)` and `ticker.standard(symbol)` as in the reference examples.
- Type system: "Comparison operators always return false if at least one operand is na", and using `na` directly as an operand of `==` is a compilation error; the script tests with `na()`. Migration guide to v6: "bool values can no longer be na", and "the and and or conditions are now evaluated lazily", with the fix "extract calls that rely on historical context to the global scope", which is why `ta.atr` and `ta.stdev` are called at the global scope on every bar and only their values are read inside the conditions.

## License

MIT. See `LICENSE`.
