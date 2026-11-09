"""Frozen B100 Type-7 CI inversion extended to B1000, and existing Fieller.

No models, random draws, participants, network or I/O. The alpha grid remains
j / actual_valid_B, strict null exclusion; this is not an exact tail P.
"""
import math
import statistics

METHOD = 'TWO_SIDED_FROZEN_TYPE7_PERCENTILE_CI_INVERSION_ALPHA_GRID_1_OVER_VALID_B'
VERSION = 'B1000_16_20261028_v1'
ALPHA = .05
Z95 = 1.959963984540054


def number(value):
    if value is None or isinstance(value, bool):
        return None
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def q7(values, probability):
    x = sorted(values)
    if not x or not 0 <= probability <= 1 or not all(math.isfinite(z) for z in x):
        raise ValueError('FINITE_NONEMPTY_SAVED_DISTRIBUTION_REQUIRED')
    h = (len(x) - 1) * probability
    j = int(h)
    g = h - j
    return x[j] + g * (x[min(j + 1, len(x) - 1)] - x[j])


def invert(values, point):
    """Legacy search and interpolation unchanged; actual-B caveat replaces B100."""
    if len(values) < 2 or not all(math.isfinite(v) for v in values) or number(point) is None:
        raise ValueError('INSUFFICIENT_OR_NONFINITE_ESTIMATES')
    n = len(values)
    low, high = q7(values, .025), q7(values, .975)
    p, previous, found = 1., (n - 1) / n, False
    # Sorting once gives exactly the same Type-7 arithmetic as the frozen
    # repeated-sort search, without changing the alpha grid or stopping rule.
    ordered = sorted(values)
    def quantile(probability):
        h = (n - 1) * probability
        j = int(h)
        return ordered[j] + (h-j)*(ordered[min(j+1,n-1)]-ordered[j])
    for j in range(1, n + 1):
        a = j / n
        if quantile(a / 2) > 0 or quantile(1 - a / 2) < 0:
            p, previous, found = a, (j - 1) / n, True
            break
    if not found:
        previous = None
    mean = sum(values) / n
    sd = math.sqrt(sum((v - mean) ** 2 for v in values) / (n - 1))
    return dict(P_raw=p, P_method=METHOD, P_method_version=VERSION,
        P_value_type='UPWARD_GRID_CI_INVERSION_SCREENING_NOT_EXACT_TAIL_PROBABILITY',
        P_grid_step=1/n, P_grid_first_step=p == 1/n,
        P_search_lower=previous, P_search_upper=p,
        P_search_bracket='[lower,upper);EMPIRICAL_CI_INVERSION_INFIMUM_ONLY' if found else 'NO_STRICT_EXCLUSION_ON_GRID;P_SET_TO_1',
        P_tail_unobserved=min(values)>0 or max(values)<0,
        P_tail_MC_upper95_if_iid=(1-.05**(1/n)) if min(values)>0 or max(values)<0 else None,
        P_MC_caveat='SEARCH_BRACKET_IS_NOT_TRUE_P_BOUND;FINITE_ACTUAL_B_TAIL_PRECISION_LIMITED;NO_POPULATION_P_ZERO_CLAIM',
        CI_excludes_zero=low>0 or high<0, CI95_type7_low=low, CI95_type7_high=high,
        P_reject_at_05=p<=ALPHA, CI_P_05_agrees=((p<=ALPHA)==(low>0 or high<0)),
        bootstrap_mean=mean, bootstrap_SD=sd, bootstrap_bias=mean-point,
        bootstrap_bias_over_SD=(mean-point)/sd if sd else None,
        P_bias_correction='NONE_MATCHES_FROZEN_PERCENTILE;NOT_BCA_NOT_STUDENTIZED',
        tested_B=n, theta_null=0, alternative='two.sided')


def bh(values, family_n):
    """R BH n=m_planned, actual missing raw P stays missing, ties order-invariant."""
    clean = [number(p) for p in values]
    observed = sorted((p, i) for i,p in enumerate(clean) if p is not None)
    if type(family_n) is not int or family_n < len(observed) or any(not 0<p<=1 for p,_ in observed):
        raise ValueError('INVALID_BH_FAMILY_OR_P')
    out, running = [None]*len(values), 1.
    for rank in range(len(observed),0,-1):
        p,i=observed[rank-1]
        running=min(running,family_n*p/rank)
        out[i]=running
    return out


def paired_covariance(numerator, denominator):
    if len(numerator)!=len(denominator) or len(numerator)<2:
        raise ValueError('PAIRED_ESTIMATES_REQUIRED')
    if any(number(x) is None for x in list(numerator)+list(denominator)):
        raise ValueError('FINITE_PAIRED_ESTIMATES_REQUIRED')
    xbar,ybar=statistics.fmean(numerator),statistics.fmean(denominator)
    dx,dy=[x-xbar for x in numerator],[y-ybar for y in denominator]
    scale=len(dx)-1
    return (math.fsum(x*x for x in dx)/scale, math.fsum(y*y for y in dy)/scale,
            math.fsum(x*y for x,y in zip(dx,dy))/scale)


def fieller_set(ie, te, var_ie, var_te, covariance, z=Z95):
    """Existing solver, sample covariance NOT divided by B; unbounded retained."""
    if any(number(x) is None for x in (ie,te,var_ie,var_te,covariance,z)) or z<=0:
        raise ValueError('FINITE_FIELLER_INPUT_REQUIRED')
    if var_ie<0 or var_te<0:
        raise ValueError('NONNEGATIVE_VARIANCES_REQUIRED')
    product=var_ie*var_te
    if covariance*covariance-product>1e-12*max(product,covariance*covariance,1e-300):
        raise ValueError('POSITIVE_SEMIDEFINITE_COVARIANCE_REQUIRED')
    a=te*te-z*z*var_te
    b=-2*(ie*te-z*z*covariance)
    c=ie*ie-z*z*var_ie
    scale=max(abs(a),abs(b),abs(c))
    aa,bb,cc=(a/scale,b/scale,c/scale) if scale else (0.,0.,0.)
    intervals=[]
    if aa==0:
        if bb==0:
            kind='ALL_REAL' if cc<=0 else 'EMPTY'
            intervals=[[-math.inf,math.inf]] if cc<=0 else []
        else:
            root=-cc/bb
            kind='HALF_LINE_LEFT' if bb>0 else 'HALF_LINE_RIGHT'
            intervals=[[-math.inf,root]] if bb>0 else [[root,math.inf]]
    else:
        discriminant=bb*bb-4*aa*cc
        tolerance=32*math.ulp(max(abs(bb*bb),abs(4*aa*cc),1e-300))
        if -tolerance<=discriminant<0:
            discriminant=0.
        if discriminant<0:
            kind='EMPTY' if aa>0 else 'ALL_REAL'
            intervals=[] if aa>0 else [[-math.inf,math.inf]]
        elif discriminant==0:
            root=-bb/(2*aa)
            kind='BOUNDED' if aa>0 else 'ALL_REAL'
            intervals=[[root,root]] if aa>0 else [[-math.inf,math.inf]]
        else:
            q=-.5*(bb+math.copysign(math.sqrt(discriminant),bb))
            roots=sorted((q/aa,cc/q))
            kind='BOUNDED' if aa>0 else 'DISJOINT_UNBOUNDED'
            intervals=[roots] if aa>0 else [[-math.inf,roots[0]],[roots[1],math.inf]]
    def token(x):
        return '-Inf' if x==-math.inf else 'Inf' if x==math.inf else x
    def label(x):
        return str(token(x)) if not math.isfinite(x) else format(x,'.12g')
    text='empty' if not intervals else ' U '.join(
        ('(' if lo==-math.inf else '[')+label(lo)+', '+label(hi)+(')' if hi==math.inf else ']')
        for lo,hi in intervals)
    return dict(kind=kind,intervals=[[token(lo),token(hi)] for lo,hi in intervals],text=text,
        A=a,B=b,C=c,critical=z,covariance_definition='SAMPLE_COVARIANCE_OF_PAIRED_BOOTSTRAP_ESTIMATES_NOT_DIVIDED_BY_B')

