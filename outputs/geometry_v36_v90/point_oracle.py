"""V13 adaptive full-operator bounds; V14 exact reflection/shift orbit reuse."""
from fractions import Fraction as F
import global_search as g


def width(record):
    return F(record['rayleigh_quotient'])-F(record['spectral_certificate']['lower'])


class Oracle:
    def __init__(self,p,tolerance,start_modes=4,max_modes=24,adaptive=True,orbits=False):
        if not 4<=start_modes<=max_modes<=32:raise ValueError('Require 4 <= start <= max <= 32')
        if not F(1,10**30)<=tolerance<=1:raise ValueError('Invalid point tolerance')
        self.p=p;self.tolerance=F(tolerance);self.start=start_modes;self.maximum=max_modes
        self.adaptive=adaptive;self.orbits=orbits;self.records=[];self.cache={};self.orbit_cache={}
        self.solves=0;self.reuses=0;self.attempts=[]

    def _compute(self,q,tolerance,start):
        chosen=None
        sizes=list(range(start,self.maximum+1,2))
        if sizes[-1]!=self.maximum:sizes.append(self.maximum)
        for n in sizes:
            sc=g.v4.guided(q,1,n,min(F(1),tolerance/2),range_proof=g.v5.make_range(q))['certificate']
            c=g.v3.old.polynomial_trial(q,n,even=not any(q[1::2]))
            mu=g.v3.old.residual_statistics(q,c)[0]
            chosen={'spectral_certificate':sc,'trial_legendre':list(map(str,c)),'rayleigh_quotient':str(mu)}
            self.solves+=1;self.attempts.append({'modes':n,'width':str(width(chosen)),'target':str(tolerance)})
            if width(chosen)<=tolerance or not self.adaptive:break
        return chosen

    @staticmethod
    def orbit(q):
        direct=tuple(q[1:]);reflected=tuple((-v if i%2==0 else v) for i,v in enumerate(q[1:]))
        flip=reflected<direct
        return (F(0),)+(reflected if flip else direct),flip,q[0]

    def _record(self,x,tolerance,start):
        q=g.at(self.p,x)
        if not self.orbits:return {'at':list(map(str,x)),**self._compute(q,tolerance,start)}
        key,flip,shift=self.orbit(q)
        canonical=self.orbit_cache.get(key)
        if canonical is None or (self.adaptive and width(canonical)>tolerance):
            canonical=self._compute(list(key),tolerance,start);self.orbit_cache[key]=canonical
        else:self.reuses+=1
        sc=canonical['spectral_certificate']
        kernel=g.v3.Kernel(q,sc['modes'],sc['tail_policy'],g.v5.make_range(q))
        transformed=g.v3.certificate(kernel,1,F(sc['lower'])+shift,F(sc['upper'])+shift,sc['upper_kind'])
        c=[F(z)*(-1 if flip and j%2 else 1) for j,z in enumerate(canonical['trial_legendre'])]
        mu=g.v3.old.residual_statistics(q,c)[0]
        return {'at':list(map(str,x)),'spectral_certificate':transformed,'trial_legendre':list(map(str,c)),'rayleigh_quotient':str(mu)}

    def point(self,x):
        x=tuple(map(F,x))
        if x in self.cache:return self.cache[x]
        i=len(self.records);self.records.append(self._record(x,self.tolerance,self.start));self.cache[x]=i
        return i

    def refine(self,i,target):
        old=self.records[i];n=old['spectral_certificate']['modes']
        if width(old)<=target:return False
        start=min(self.maximum,n+2)
        candidate=self._record(tuple(map(F,old['at'])),F(target),start)
        # Keep the better lower endpoint and trial independently at the same point.
        if F(candidate['spectral_certificate']['lower'])<F(old['spectral_certificate']['lower']):
            candidate['spectral_certificate']=old['spectral_certificate']
        if F(candidate['rayleigh_quotient'])>F(old['rayleigh_quotient']):
            candidate['trial_legendre']=old['trial_legendre'];candidate['rayleigh_quotient']=old['rayleigh_quotient']
        self.records[i]=candidate
        return width(candidate)<width(old)

    def statistics(self):
        return {'points':len(self.records),'spectral_searches':self.solves,'orbit_reuses':self.reuses,
                'modes_histogram':{str(n):sum(r['spectral_certificate']['modes']==n for r in self.records)
                                   for n in sorted({r['spectral_certificate']['modes'] for r in self.records})},
                'maximum_point_width':str(max(map(width,self.records)))}
