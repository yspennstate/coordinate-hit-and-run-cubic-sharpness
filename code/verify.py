"""Single-thread checks for the rotated box bipyramid in main.tex."""
import json, math
from fractions import Fraction as F
from pathlib import Path
import numpy as np


def integral_polynomial(d, coeff):
    return sum((F(c, d+k) for k,c in enumerate(coeff)), F(0))


def exact_checks():
    checks=[]
    for n in range(2,101):
        d=n-2
        moment=F(1) if d==0 else integral_polynomial(d,[1,-2,0,2,-1])/integral_polynomial(d,[1,0,-1])
        expected=F(6*n,(n-1)*(n+1)*(n+2))
        assert moment==expected,(n,moment,expected)
        # Independent y-slice integration: t is uniform in [-L*r,L*r].
        vt=F((2*n)**2,3)*F(n,n+2)
        rayleigh=F(2*4,3*n)*moment/vt
        assert rayleigh==F(12,n**3*(n*n-1))
        assert F(n*n+n-1,n*n-1)<=F(5,3)
        # Each facet normal has norm at most 1; hence B_2 is contained.
        assert F(1,4)+F(1,4)<=1
        assert F(1,4)+F(1,4*n*n)<=1
        checks.append(n)
    return {'dimensions':checks,'rational_moment_and_bound_checks':'passed'}


def chord(y,t,m,a,L,s):
    # Fixing x2 gives t=y+c (s=1), fixing x1 gives t=-y+c (s=-1).
    c=t-s*y
    lo=-a*(1-m); hi=a*(1-m)
    for ey in [-1,1]:
        for et in [-1,1]:
            b=ey/a+et*s/L
            rhs=1-et*c/L
            if b>0:hi=np.minimum(hi,rhs/b)
            elif b<0:lo=np.maximum(lo,rhs/b)
            else:assert np.min(rhs)>=-1e-10
    assert np.max(lo-y)<1e-10
    assert np.max(y-hi)<1e-10
    assert np.min(hi-lo)>-1e-10
    for yy in [lo,hi,(lo+hi)/2]:
        assert np.max(np.abs(yy)/a+np.abs(s*yy+c)/L)<=1+1e-10
        assert np.max(np.abs(yy)/a+m)<=1+1e-10
    return (hi-lo)**2/12


def numerical(n,N=30000):
    rng=np.random.default_rng(410000+n)
    a=R=2.;L=2.*n
    r=rng.random(N)**(1/n)
    y=a*(1-r)*rng.choice([-1,1],N)
    t=L*r*rng.uniform(-1,1,N)
    if n>2:
        z=R*r[:,None]*rng.uniform(-1,1,(N,n-2))
        m=np.max(np.abs(z),axis=1)/R
    else:m=np.zeros(N)
    assert np.max(np.abs(y)/a+np.maximum(np.abs(t)/L,m))<1+1e-12
    v1=chord(y,t,m,a,L,1)
    v2=chord(y,t,m,a,L,-1)
    envelope=a*a*(1-m)**2/3
    ratio=np.maximum(v1,v2)/envelope
    assert np.max(ratio)<=1+1e-9
    exact_m=6*n/((n-1)*(n+1)*(n+2))
    exact_v=L*L*n/(3*(n+2))
    measured_D=(v1+v2)/n
    measured_Q=np.mean(measured_D)/exact_v
    bound_Q=12/(n**3*(n*n-1))
    return {'n':n,'samples':N,'seed':410000+n,
      'moment_exact':exact_m,'moment_measured':float(np.mean((1-m)**2)),
      'moment_standard_error':float(np.std((1-m)**2,ddof=1)/math.sqrt(N)),
      'variance_t_exact':exact_v,'second_moment_t_measured':float(np.mean(t*t)),
      'rayleigh_measured':float(measured_Q),
      'rayleigh_standard_error':float(np.std(measured_D,ddof=1)/math.sqrt(N)/exact_v),
      'rayleigh_upper_bound':bound_Q,
      'rayleigh_over_bound':float(measured_Q/bound_Q),
      'maximum_chord_variance_over_envelope':float(np.max(ratio)),
      'n3_times_gap_CPI_upper_bound':192*(n*n+n-1)/(math.pi**2*(n*n-1))}


def boundary_checks():
    # Feasible vertices, base corners and central slices, including n=2's empty z.
    count=0
    for n in range(2,9):
        a=2.;L=2.*n
        for y in [-a,-a+1e-10,-1.,0.,1.,a-1e-10,a]:
            r=max(0.,1-abs(y)/a)
            for t in [-L*r,0.,L*r]:
                for m in ([0.] if n==2 else [0.,r/2,r]):
                    for s in [-1,1]:
                        v=chord(np.array([y]),np.array([t]),np.array([m]),a,L,s)[0]
                        assert v<=a*a*(1-m)**2/3+1e-10
                        count+=1
    return {'tested_boundary_slices':count,'dimensions':list(range(2,9)),'result':'passed'}


def controls():
    # Without rotation, resampling t has conditional mean zero and all other
    # updates retain t. Its eigenvalue is 1-1/n; the full gap is not identified.
    # Without taper (rotated rectangle times cube), the diagonal strip
    # has a chord scale independent of n; test the interior fiber exactly.
    out=[]
    for n in [4,16,64]:
        a=2.;L=2.*n
        # Untapered rectangle: |y|<=a and |t|<=L, with t=y on this fiber.
        lo=max(-a,-L);hi=min(a,L)
        v=(hi-lo)**2/12
        assert abs(v-a*a/3)<1e-12
        out.append({'n':n,'untapered_central_conditional_variance_t':float(v),
                    'unrotated_t_eigenvalue_gap':1/n})
    return out

if __name__=='__main__':
    result={'exact':exact_checks(),'boundary':boundary_checks(),'controls':controls(),
      'numerical':[numerical(n) for n in [2,3,4,8,16,32,64,128]],
      'warning':'Monte Carlo estimates are Rayleigh quotients, not spectral gaps.'}
    Path('verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Exact rational dimensions 2..100 and boundary slices passed.')
    print(' n     Rayleigh/bound    n^3*(gap*CPI) bound')
    for r in result['numerical']:
        print('%3d    %.8f           %.8f'%(r['n'],r['rayleigh_over_bound'],r['n3_times_gap_CPI_upper_bound']))
