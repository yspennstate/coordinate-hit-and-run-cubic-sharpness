"""Rational chord integration adapted from the independent construction check."""
from fractions import Fraction as F
from pathlib import Path
import json

def vertices(m, a, L):
    h=1-m
    return [(0,L),(a*h,L*m),(a*h,-L*m),(0,-L),(-a*h,-L*m),(-a*h,L*m)]

def area(v):
    return abs(sum(y*t1-t*y1 for (y,t),(y1,t1) in zip(v,v[1:]+v[:1])))/2

def width(v,c):
    # Intersect the projected polygon edges with the vertical line t-y=c.
    intersections=[]
    for (y,t),(yy,tt) in zip(v,v[1:]+v[:1]):
        c0,c1=t-y,tt-yy
        if min(c0,c1)<=c<=max(c0,c1):
            if c0==c1:
                intersections.extend([F(y),F(yy)])
            else:
                intersections.append(y+(yy-y)*(c-c0)/(c1-c0))
    assert intersections
    return max(intersections)-min(intersections)

def slice_integrals(m,a,L):
    v=vertices(m,a,L)
    cuts=sorted(set(t-y for y,t in v))
    i1=i3=F(0)
    for lo,hi in zip(cuts,cuts[1:]):
        # Inside each interval the width is affine. Evaluate away from vertices.
        w1=width(v,(2*lo+hi)/3); w2=width(v,(lo+2*hi)/3)
        w0=2*w1-w2; w3=2*w2-w1
        assert min(w0,w3)>=0
        i1+=(hi-lo)*(w0+w3)/2
        i3+=(hi-lo)*(w0**3+w0*w0*w3+w0*w3*w3+w3**3)/4
    assert i1==area(v)==2*a*L*(1-m*m)
    return i1,i3

def multiply(p,q):
    out=[F(0)]*(len(p)+len(q)-1)
    for i,x in enumerate(p):
        for j,y in enumerate(q):out[i+j]+=x*y
    return out

def interpolate(xs,ys):
    out=[F(0)]*len(xs)
    for i,(x,y) in enumerate(zip(xs,ys)):
        basis=[F(1)]; den=F(1)
        for j,z in enumerate(xs):
            if i!=j:basis=multiply(basis,[-z,F(1)]);den*=x-z
        for j,b in enumerate(basis):out[j]+=y*b/den
    return out

def evaluate(p,x):
    return sum((a*x**i for i,a in enumerate(p)),F(0))

def integrate_weighted(p,lo,hi,d):
    return sum((a*(hi**(d+i)-lo**(d+i))/(d+i) for i,a in enumerate(p)),F(0))

def exact_rayleigh(n,a=F(2),L=None):
    L=F(2*n) if L is None else L
    d=n-2
    if not d:
        A,I=slice_integrals(F(0),a,L)
        ev=I/(12*A)
    else:
        # All projected vertex-order changes, including L<a and L=a.
        cs0=[t-y for y,t in vertices(F(0),a,L)]
        cs1=[t-y for y,t in vertices(F(1),a,L)]
        lines=[(b,e-b) for b,e in zip(cs0,cs1)]
        breaks={F(0),F(1)}
        for b,s in lines:
            for bb,ss in lines:
                if s!=ss:
                    root=(bb-b)/(s-ss)
                    if 0<root<1:breaks.add(root)
        breaks=sorted(breaks)
        total=F(0)
        for lo,hi in zip(breaks,breaks[1:]):
            xs=[lo+(hi-lo)*F(i,6) for i in range(1,6)]
            p=interpolate(xs,[slice_integrals(x,a,L)[1] for x in xs])
            # Validate quartic interpolation at two extra rational points.
            for u in [F(1,7),F(6,7)]:
                x=lo+(hi-lo)*u
                assert evaluate(p,x)==slice_integrals(x,a,L)[1]
            total+=integrate_weighted(p,lo,hi,d)
        ev=total/(12*(2*a*L)*(F(1,d)-F(1,d+2)))
    # Order-statistic representation: rho ~ Beta(n-1,2), maximum index uniform.
    erho2=F((n-1)*n,(n+1)*(n+2))
    vt=L*L*erho2*(F(1,n-1)+F(n-2,3*(n-1)))
    assert vt==L*L*F(n,3*(n+2))
    # Derive moments via a different disintegration over box shells.
    if d:
        erho1=F(n-1,n+1)
        em=erho1*(F(d,n-1)+F(1,n-1)*F(d,d+1))
        em2=erho2*(F(d,n-1)+F(1,n-1)*F(d,d+2))
        slack=1-2*em+em2
    else:slack=F(1)
    assert slack==F(6*n,(n-1)*(n+1)*(n+2))
    quotient=F(2,n)*ev/vt
    bound=12*a*a/(L*L*n*(n*n-1))
    assert 0<quotient<=bound
    return {'n':n,'a':str(a),'L':str(L),'Var_t':str(vt),
            'E_slack_squared':str(slack),'E_conditional_variance':str(ev),
            'exact_Rayleigh':str(quotient),'exact_Rayleigh_float':float(quotient),
            'envelope':str(bound),'ratio':float(quotient/bound)}

def boundary_cases():
    count=0
    for a,L in [(F(2),F(4)),(F(2),F(2)),(F(1,100),F(10)),(F(3),F(6))]:
        for m in [F(0),F(1,1000000),F(1,2),F(999999,1000000),F(1)]:
            v=vertices(m,a,L)
            cuts=sorted(set(t-y for y,t in v))
            test=cuts+[(u+v)/2 for u,v in zip(cuts,cuts[1:])]
            for c in test:
                ell=width(v,c)
                assert 0<=ell<=2*a*(1-m)
                count+=1
    return count

if __name__ == '__main__':
    rows=[exact_rayleigh(n) for n in [2,3,4,8,16,32,64,128]]
    generic=[exact_rayleigh(n,a,L) for n in [2,3,5,10]
             for a,L in [(F(3),F(6)),(F(2),F(2)),(F(6),F(3)),
                         (F(1,100),F(10)),(F(10),F(1,100))]]
    boundaries=boundary_cases()
    assert rows[0]['exact_Rayleigh']=='5/27'
    assert rows[1]['exact_Rayleigh']=='17/512'
    result={'family':rows,'varied_scales':generic,'boundary_fibers':boundaries}
    Path('exact_chords.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Exact Rayleigh quotients: n=2:',rows[0]['exact_Rayleigh'],
          '; n=3:',rows[1]['exact_Rayleigh'])
    print(len(rows),'family cases;',len(generic),'varied-scale cases;',boundaries,'boundary fibers: passed')
