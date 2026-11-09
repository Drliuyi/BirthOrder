#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <math.h>
#include <float.h>
#include <fenv.h>

/* DOUBLE_ONLY_BINARY64 realization of the approved competing-hazard
 * exponential, NOT discrete AJ. The formula and integration method are
 * unchanged; only working precision and stable evaluation order change.
 * No probability clipping, row renormalization, or survival tail cutoff.
 * The complete input is checked before the first transition. */
#if defined(__FAST_MATH__) || \
    (defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__ > 0)
#error C39_EXP_DOUBLE_ONLY_UNSAFE_FLOATING_POINT_COMPILER_MODE
#endif
#if FLT_RADIX != 2 || DBL_MANT_DIG != 53 || DBL_MAX_EXP != 1024 || \
    DBL_MIN_EXP != -1021 || FLT_EVAL_METHOD != 0
#error C39_EXP_DOUBLE_ONLY_BINARY64_EVALUATION_REQUIRED
#endif
#if defined(DBL_HAS_SUBNORM) && DBL_HAS_SUBNORM != 1
#error C39_EXP_DOUBLE_ONLY_SUBNORMAL_REPRESENTATION_REQUIRED
#endif

static void environment_guard(void) {
  if (sizeof(double) != 8 || fegetround() != FE_TONEAREST)
    error("C39_EXP_DOUBLE_ONLY_BINARY64_ROUND_TO_NEAREST_REQUIRED");
  volatile double tiny = nextafter(0.0, 1.0), one = 1.0, two = 2.0;
  volatile double smallest_normal = DBL_MIN;
  if (!(tiny > 0.0 && tiny * one == tiny && tiny + tiny > tiny &&
        smallest_normal / two > 0.0 && (tiny + tiny) / two == tiny))
    error("C39_EXP_DOUBLE_ONLY_GRADUAL_UNDERFLOW_REQUIRED");
}

/* Existing exported symbol/API retained; reports the working mantissa. */
SEXP c39_exp_environment(void) {
  environment_guard();
  return ScalarInteger(DBL_MANT_DIG);
}

static void add_compensated(double value, double *sum, double *err) {
  double y = value - *err, t = *sum + y;
  *err = (t - *sum) - y;
  *sum = t;
}

typedef struct {
  double mantissa;
  int exponent;
} scaled_positive;

/* The product is kept as mantissa * 2^exponent without overflow or
 * premature subnormal multiplication. All arithmetic remains binary64. */
static scaled_positive scaled_product(double h, double rr) {
  scaled_positive p;
  if (h == 0.0) {
    p.mantissa = 0.0;
    p.exponent = 0;
  } else {
    int eh, er;
    p.mantissa = frexp(h, &eh) * frexp(rr, &er);
    p.exponent = eh + er;
    if (p.mantissa < 0.5) {
      p.mantissa *= 2.0;
      p.exponent--;
    }
  }
  return p;
}

/* A conservative bound leaves an entire exponent of headroom even if the
 * division used to check the bound rounds upward. It selects arithmetic,
 * never a probability or survival approximation. */
static int moderate_product(double h, double rr) {
  double bound = DBL_MAX * 0.25;
  return rr <= 1.0 ? h <= bound : h <= bound / rr;
}

/* Return 1 only when q can be represented; never form an overflowing q.
 * The common ordinary-input path retains direct product/sum evaluation. */
static int finite_transition(double h1, double h2, double rr1, double rr2,
                             double *q, double *w1, double *w2) {
  if (moderate_product(h1, rr1) && moderate_product(h2, rr2)) {
    double a = h1 * rr1, b = h2 * rr2;
    *q = a + b;
    if (*q > 0.0) {
      *w1 = a / *q;
      *w2 = b / *q;
    }
    return 1;
  }
  scaled_positive a = scaled_product(h1, rr1);
  scaled_positive b = scaled_product(h2, rr2);
  if (a.mantissa == 0.0 && b.mantissa == 0.0) {
    *q = 0.0;
    return 1;
  }
  int exponent = a.mantissa == 0.0 ? b.exponent :
                 b.mantissa == 0.0 ? a.exponent :
                 a.exponent > b.exponent ? a.exponent : b.exponent;
  double x = a.mantissa == 0.0 ? 0.0 :
             scalbn(a.mantissa, a.exponent - exponent);
  double y = b.mantissa == 0.0 ? 0.0 :
             scalbn(b.mantissa, b.exponent - exponent);
  double total = x + y;
  if (total >= 1.0) {
    total *= 0.5;
    exponent++;
  }
  /* 0x1.fffffffffffffp-1 is DBL_MAX scaled by 2^-DBL_MAX_EXP. */
  if (exponent > DBL_MAX_EXP ||
      (exponent == DBL_MAX_EXP && total > 0x1.fffffffffffffp-1))
    return 0;
  *q = scalbn(total, exponent);
  if (*q > 0.0) {
    *w1 = a.mantissa == 0.0 ? 0.0 :
          scalbn(a.mantissa / total, a.exponent - exponent);
    *w2 = b.mantissa == 0.0 ? 0.0 :
          scalbn(b.mantissa / total, b.exponent - exponent);
  }
  return 1;
}

static double survival_step(double s, double q, double remain) {
  if (remain > 0.0 && remain < DBL_MIN) {
    /* Same s * exp(-q), split before the subnormal rounding of exp(-q).
     * Example: a rounded minimum-subnormal remain must not be multiplied
     * as though that rounded value were the exact exponential. No cutoff. */
    double half = exp(-0.5 * q);
    return (s * half) * half;
  }
  return s * remain;
}

SEXP c39_piecewise_exponential(SEXP H, SEXP RR) {
  environment_guard();
  if (TYPEOF(H) != REALSXP || TYPEOF(RR) != REALSXP ||
      !isMatrix(H) || !isMatrix(RR))
    error("C39_EXP_DOUBLE_MATRICES_REQUIRED");
  SEXP hd = getAttrib(H, R_DimSymbol), rd = getAttrib(RR, R_DimSymbol);
  int m = INTEGER(hd)[0], n = INTEGER(rd)[0];
  if (INTEGER(hd)[1] != 2 || INTEGER(rd)[1] != 2 || n < 1)
    error("C39_EXP_DIMENSION_INVALID");
  const double *h = REAL(H), *rr = REAL(RR);
  R_xlen_t hm = (R_xlen_t)m, rn = (R_xlen_t)n;
  for (R_xlen_t j = 0; j < XLENGTH(H); j++)
    if (!R_FINITE(h[j]) || h[j] < 0.0)
      error("C39_EXP_INVALID_BASELINE index=%lld", (long long)j + 1);
  for (R_xlen_t i = 0; i < XLENGTH(RR); i++)
    if (!R_FINITE(rr[i]) || rr[i] <= 0.0)
      error("C39_EXP_INVALID_RR index=%lld", (long long)i + 1);
  /* Exact exponential semigroup compression, not time-grid approximation:
   * consecutive jumps of one and the SAME cause commute at fixed individual
   * RR. Mixed-cause jumps and cause changes are immutable barriers. Validate
   * all original inputs above, and retain them for the overflow fallback.
   * A sum that cannot be represented is never published as an Inf baseline. */
  double *ch = (double *) R_alloc((size_t)(m ? 2*m : 1), sizeof(double));
  int cm = 0;
  for (int j=0; j<m;) {
    double a=h[j], b=h[j+hm];
    int cause = a>0.0 && b==0.0 ? 1 : b>0.0 && a==0.0 ? 2 : 0;
    int k=j+1;
    if (cause) {
      double sum = cause==1 ? a : b, err=0.0;
      while (k<m && (cause==1 ? h[k]>0.0 && h[k+hm]==0.0 : h[k+hm]>0.0 && h[k]==0.0)) {
        double v=cause==1 ? h[k] : h[k+hm];
        double y=v-err;
        if (y>0.0 && sum>DBL_MAX-y) break;
        double total=sum+y;
        if (!R_FINITE(total)) break;
        err=(total-sum)-y; sum=total; k++;
      }
      if (cause==1) a=sum; else b=sum;
    }
    ch[cm]=a; ch[cm+m]=b; cm++; j=k;
  }
  double max1=0.0, max2=0.0;
  for (int j=0;j<cm;j++) {max1=fmax(max1,ch[j]);max2=fmax(max2,ch[j+m]);}
  /* p(3), mass residual, incidence/survival relation residual,
   * log-scaled overflow steps, underflowed survival indicator. */
  SEXP out = PROTECT(allocMatrix(REALSXP, n, 7));
  double *o = REAL(out);
  for (int i = 0; i < n; i++) {
    if (i % 32 == 0) R_CheckUserInterrupt();
    double s = 1.0, f1 = 0.0, f2 = 0.0, e1 = 0.0, e2 = 0.0;
    double scaled = 0.0;
    int underflow = 0;
    /* The fast row is proven free of overflow, including merged jumps.
     * Otherwise use ALL original transitions, preserving the original count
     * and meaning of the overflow diagnostic (output column 6). */
    int fast = moderate_product(max1,rr[i]) && moderate_product(max2,rr[i+rn]);
    int steps=fast ? cm : m;
    const double *hh=fast ? ch : h;
    for (int j = 0; j < steps; j++) {
      double q = 0.0, w1 = 0.0, w2 = 0.0, event, remain;
      int representable;
      double h1=hh[j],h2=hh[j+(fast ? m : hm)];
      if (fast) {
        double a=h1*rr[i],b=h2*rr[i+rn]; q=a+b;representable=1;
        if(q>0.0) { w1=a/q;w2=b/q; }
      } else representable=finite_transition(h1,h2,rr[i],rr[i+rn],&q,&w1,&w2);
      if (representable) {
        if (q == 0.0) continue; /* Exact zero or representational underflow. */
        /* Compute the SMALLER tail directly; subtract only the larger one.
         * Small-q incidence retains expm1 accuracy, large-q survival retains
         * exp accuracy. No cutoff or clipping; survival_step is unchanged. */
        if(q<=0x1.62e42fefa39efp-1) {
          event=-expm1(-q);remain=1.0-event;
        } else { remain=exp(-q);event=1.0-remain; }
      } else {
        /* Finite inputs can have an unrepresentable total hazard. Its q
         * is never formed. LSE preserves the competing-intensity ratio.
         * event=1, remain=0 are the binary64 representations of this
         * exponential, NOT a survival-tail rule. */
        double la = h1 > 0.0 ? log(h1) + log(rr[i]) : -INFINITY;
        double lb = h2 > 0.0 ? log(h2) + log(rr[i + rn]) : -INFINITY;
        double mx = fmax(la, lb), x = exp(la - mx), y = exp(lb - mx);
        w1 = x / (x + y);
        w2 = y / (x + y);
        event = 1.0;
        remain = 0.0;
        scaled++;
      }
      if (!(isfinite(w1) && isfinite(w2) && isfinite(event) && isfinite(remain)) ||
          w1 < 0.0 || w1 > 1.0 || w2 < 0.0 || w2 > 1.0 ||
          event < 0.0 || event > 1.0 ||
          remain < 0.0 || remain > 1.0)
        error("C39_EXP_NONFINITE_OR_INVALID_TRANSITION row=%d jump=%d", i + 1, j + 1);
      add_compensated(s * event * w1, &f1, &e1);
      add_compensated(s * event * w2, &f2, &e2);
      double before = s;
      s = representable ? survival_step(s, q, remain) : s * remain;
      if (before > 0.0 && s == 0.0) underflow = 1;
    }
    double p1 = f1, p2 = f2, ps = s;
    if (!R_FINITE(p1) || !R_FINITE(p2) || !R_FINITE(ps) ||
        p1 < 0.0 || p1 > 1.0 || p2 < 0.0 || p2 > 1.0 || ps < 0.0 || ps > 1.0)
      error("C39_EXP_INVALID_PROBABILITY_NO_CLIP row=%d", i + 1);
    o[i] = p1;
    o[i + rn] = p2;
    o[i + 2 * rn] = ps;
    o[i + 3 * rn] = fabs((p1 + p2 + ps) - 1.0);
    o[i + 4 * rn] = fabs((p1 + p2) - (1.0 - ps));
    o[i + 5 * rn] = scaled;
    o[i + 6 * rn] = underflow;
  }
  UNPROTECT(1);
  return out;
}
