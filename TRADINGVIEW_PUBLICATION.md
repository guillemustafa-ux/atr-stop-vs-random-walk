# ATR Stops vs a Random Walk: How Often Noise Hits Them

Description to paste as is (between the rules):

---
WHAT IT DOES
A stop placed a few ATRs from price is meant to sit outside the noise. This indicator counts how often price reaches it anyway, and compares that count with what a pure random walk of the same volatility would give. At regular entries that never overlap, it places a long stop and a short stop at 1, 1.5, 2 and 3 ATRs from the close and checks the next N bars: the lows for the long stop, the highs for the short stop. For each multiple the table on the last bar shows the long, short and pooled hit rates with a 95% confidence interval, the random-walk rate for the same stops and the same volatility, the difference in points and a verdict.

WHAT IS ORIGINAL HERE
ATR stop and trailing-stop tools draw levels; they do not say how often price reaches them, or how often it would by chance alone. Here the chance has a closed form. For driftless Brownian motion in log price with volatility sigma per bar, the probability of touching a level at log distance d within N bars is 2 x Phi(-d / (sigma x sqrt(N))), by the reflection principle. Each resolved event gets that prediction from its own entry, ATR and sigma, so the measured rate and the random-walk rate come from the same bars. Pine has no normal distribution function, so Phi is built from the Abramowitz and Stegun 7.1.26 approximation of erf, with an error below 1e-7. A volatility row adds a check that does not depend on the stops: the expected range of a Brownian bar is about 1.6 sigma, so an ATR / sigma ratio far from 1.6 says the bars do not move like a random walk.

HOW IT IS MEASURED
- Entry: the close of a confirmed bar with no event open, once the ATR and the volatility window are ready. An event resolves N bars after its entry (20 by default), and the next one can start on that bar, so entries are at least N bars apart.
- ATR: ta.atr(14) on the entry bar. Sigma: the sample standard deviation of the 100 one-bar log returns that end at the entry bar.
- Stops: entry - k x ATR and entry + k x ATR; a hit is a touch by a low or a high of the N bars after the entry.
- Prediction: d = log distance from the entry to the stop, rate = 2 x Phi(-d / (sigma x sqrt(N))), averaged over both sides of all events.
- Interval: the event is the unit, because its long and short outcomes share the same bars.
- Counters move only on confirmed bars; the event still open on the last bar is not counted. No request.security call.

HOW TO USE IT
1. Keep the defaults first, then set the horizon to the bars you would hold a position.
2. Read the "difference" column: the pooled hit rate minus the random-walk rate, in points.
3. Read the verdict: "like a random walk" when the random-walk rate lies inside the interval, "hit more than a random walk" when it lies below it, "hit less than a random walk" when it lies above it; "(small sample)" under 30 events.
4. Read the volatility row: with an ATR / sigma ratio far from 1.6 the bars are not random-walk bars, and the verdicts reflect that too.
5. Set Stop 1 to 0.1: that row should read near 100%. Set the horizon to 1: every row should be small.

HOW TO READ THE TABLE
Header: symbol, timeframe, ATR length, volatility window, horizon.
One row per stop: events, long hit, short hit, pooled hit, 95% CI, random walk, difference, verdict.
Volatility: mean ATR as % of price, mean sigma per bar, their ratio, the Brownian reference 1.60, bars counted and events.
Last event: entry time, close, ATR, the lowest low and the highest high of its N bars in ATRs from the entry, and the open event that is not counted yet.

LIMITATIONS
- The null is a driftless random walk with a constant volatility per event. A trend, a change of volatility, mean reversion or fat tails can move the hit rates away from it; the verdict says that the chart differs, not why.
- Sigma comes from 100 returns and has an error of about 7% of its own, which the interval does not include.
- The ATR and sigma windows react to a change of volatility at different speeds.
- A hit is a touch of the high or the low: no fills, slippage or costs, and no order between the long and the short stop.
- Use standard candles; Heikin Ashi or Renko bars change the highs and lows.
- Entries are spaced from the first loaded bar, so history depth, symbol and timeframe change the result. Counts, not a trading rule. Not financial advice.

The chart shows BINANCE:BTCUSDT on the 1 hour timeframe with the default settings. The table is in the bottom right.
---

## Publishing checklist

- [ ] Compile `atr_stop_vs_random_walk.pine` in the Pine Editor with no errors or warnings.
- [ ] Set the chart: BINANCE:BTCUSDT, 1 hour, default inputs, table visible.
- [ ] Check the identities: no rate above 100%; the same number of events in every stop row; the pooled rate is the mean of the long and the short rates up to rounding; with the default ascending stops the long, the short and the random-walk rates never rise from one row to the next; "difference" equals pooled minus random walk up to rounding; ATR / sigma equals the ATR percentage over the sigma percentage up to rounding.
- [ ] Control: set Stop 1 (x ATR) to 0.1. The first row's pooled and random-walk rates must read near 100% (about 97% with the ratio near 1.6). Set it back to 1.
- [ ] Control: set Horizon N (bars) to 1. Every row must be small (the random walk gives about 11% at 1 ATR with the ratio near 1.6, and less further out), with about 20 times the events. Set it back to 20.
- [ ] Hand check the last event: at the circle of its entry time, the lowest low of the next 20 bars minus the entry close, divided by the ATR in the row, equals the "low" cell (0.00 when no low went below the entry).
- [ ] Check the title is at most 60 characters (it is 53).
- [ ] Publish as open source, with the MIT license text if the form allows it.
- [ ] Check that the description contains no links, contact info, or performance claims, and is between 3,000 and 4,600 characters.
- [ ] In Final touches, open Show more and fill the tags BEFORE pressing Publish; the field takes 5 tags at most, and the Tags field has to be clicked inside before each paste. Take the tags from the suggestion list and never a category name ("volatility" and "statistics" are categories, so avoid them). Suggested tags (one word each): atr, stoploss, randomwalk, brownian, hitrate.
- [ ] Before pressing Publish, add the script to the chart from the editor and read the table on screen: the snapshot shows what is on the chart at that moment. Make sure the chart has no drawings: a stray paste becomes a text drawing in the snapshot.
- [ ] Add the TradingView script URL to `README.md` after publishing.
