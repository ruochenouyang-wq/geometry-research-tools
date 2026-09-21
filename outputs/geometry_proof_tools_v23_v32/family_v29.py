"""V29: compile a polynomial Poisson equation into an all-amplitude Barta lemma."""
from fractions import Fraction as F
import hashlib,json
import error_bounds as e
import algebra as a
import v05_range as ranges


def mean(p):return sum((v/F(i+1) for i,v in enumerate(p) if i%2==0),F(0))


def laplacian(p):
    return e.add(e.mul([0,2],e.derivative(p)),e.scale(e.mul([1,0,-1],e.derivative(e.derivative(p))),-1))


def poisson(raw):
    p=e.base.potential(list(map(str,raw)));n=len(p)-1
    if not n:return [F(0)],mean(p)
    columns=[laplacian([F(0)]*j+[F(1)]) for j in range(1,n+1)]
    matrix=[[col[i] if i<len(col) else F(0) for col in columns] for i in range(1,n+1)]
    r=[F(0)]+a.solve(matrix,[-p[i] for i in range(1,n+1)])
    target=e.add(e.scale(p,-1),[mean(p)])
    if laplacian(r)!=target:raise ArithmeticError('Poisson identity failed')
    return e.trim(r),mean(p)


def reduced_polynomial(raw):
    q=e.trim(list(map(F,raw)))
    if not 1<=len(q)<=13:raise ValueError('Polynomial degree <= 12 required')
    even=not any(q[1::2]) and len(q)>1
    return (q[::2],F(0),F(1),'t_squared') if even else (q,F(-1),F(1),'t')


def range_certificate(raw,intervals,witness,epsilon):
    original=e.trim(list(map(F,raw)));q,start,end,coordinate=reduced_polynomial(original);epsilon=F(epsilon)
    if not F(1,10**12)<=epsilon<=1:raise ValueError('Range accuracy budget')
    if len(q)<=3:
        points=[start,end]
        if len(q)==3 and q[2]:
            x=-q[1]/(2*q[2])
            if start<=x<=end:points.append(x)
        location=max(points,key=lambda t:e.evaluate(q,t));upper=e.evaluate(q,location);lower=upper
        if intervals is not None or witness!=str(location):raise ValueError('Incorrect exact maximum')
    else:
        if not isinstance(intervals,list) or not 1<=len(intervals)<=256:raise ValueError('Range partition budget')
        last=start;upper=None
        for pair in intervals:
            lo,hi=map(F,pair)
            if lo!=last or not lo<hi<=end:raise ValueError('Range cover incomplete')
            val=max(ranges.interval_coefficients(q,lo,hi));upper=val if upper is None else max(upper,val);last=hi
        if last!=end:raise ValueError('Range misses endpoint')
        location=F(witness)
        if not start<=location<=end:raise ValueError('Maximum witness outside domain')
        lower=e.evaluate(q,location)
    return {'method':'polynomial_maximum','polynomial':list(map(str,original)),'coordinate':coordinate,
            'intervals':intervals,'witness':str(location),'maximum_lower':str(lower),'maximum_upper':str(upper),
            'tolerance':str(epsilon),'status':'scale_enclosed' if upper-lower<=epsilon else 'scale_gap_open'}


def maximum(q,epsilon=F(1,10**6),max_leaves=64):
    if type(max_leaves) is not int or not 1<=max_leaves<=256:raise ValueError('Range leaf budget')
    original=e.trim(list(map(F,q)));q,start,end,_=reduced_polynomial(original)
    if len(q)<=3:
        points=[start,end]
        if len(q)==3 and q[2] and start<=-q[1]/(2*q[2])<=end:points.append(-q[1]/(2*q[2]))
        return range_certificate(original,None,str(max(points,key=lambda t:e.evaluate(q,t))),epsilon)
    leaves=[(start,end,ranges.interval_coefficients(q,start,end))];samples={t:e.evaluate(q,t) for t in (start,(start+end)/2,end)}
    while max(max(v) for _,_,v in leaves)-max(samples.values())>epsilon and len(leaves)<max_leaves:
        i=max(range(len(leaves)),key=lambda i:max(leaves[i][2]));lo,hi,b=leaves.pop(i);mid=(lo+hi)/2;l,r=e.split_bernstein(b)
        leaves.extend([(lo,mid,l),(mid,hi,r)])
        for t in ((lo+mid)/2,mid,(mid+hi)/2):samples[t]=e.evaluate(q,t)
    intervals=[[str(lo),str(hi)] for lo,hi,_ in sorted(leaves)]
    return range_certificate(original,intervals,str(max(samples,key=samples.get)),epsilon)


def verify_range(cert):
    try:return cert==range_certificate(cert['polynomial'],cert['intervals'],cert['witness'],cert['tolerance'])
    except (ValueError,TypeError,KeyError,IndexError,ZeroDivisionError):return False


def assemble(directions,weight):
    dirs=[e.base.potential(list(map(str,p))) for p in directions];d=len(dirs)
    if not 1<=d<=4:raise ValueError('Family has one to four directions')
    weight=[list(map(F,row)) for row in weight]
    if len(weight)!=d or any(len(row)!=d for row in weight) or any(weight[i][j]!=weight[j][i] for i in range(d) for j in range(d)) or a.inertia(weight)!=[0,0,d]:raise ValueError('Positive definite symmetric weight required')
    inv=a.inverse(weight);solutions=[poisson(p) for p in dirs];r=[x[0] for x in solutions];means=[x[1] for x in solutions];der=[e.derivative(v) for v in r]
    g=[F(0)]
    for i in range(d):
        for j in range(d):g=e.add(g,e.scale(e.mul(der[i],der[j]),inv[i][j]))
    g=e.mul([1,0,-1],g)
    return dirs,weight,r,means,g


def build(directions,weight,range_proof,method):
    dirs,w,r,means,g=assemble(directions,weight)
    if not verify_range(range_proof) or list(map(F,range_proof['polynomial']))!=g:raise ValueError('Range is not bound to the Poisson family')
    if method=='poisson_scalar_v29' and (len(dirs)!=1 or w!=[[F(1)]]):raise ValueError('Scalar template requires one normalized weight')
    if method not in ('poisson_scalar_v29','poisson_matrix_v30'):raise ValueError('Unknown family method')
    upper=F(range_proof['maximum_upper'])
    if upper<0:raise ValueError('Nonnegative rank-one envelope required')
    cert={'method':method,'directions':[list(map(str,p)) for p in dirs],'weight':[list(map(str,row)) for row in w],
          'poisson_solutions':[list(map(str,v)) for v in r],'means':list(map(str,means)),
          'scale':str(upper),'quadratic_bound':[[str(upper*v) for v in row] for row in w],'range_proof':range_proof,
          'scope':'ground energy on the full unit sphere, all real amplitude vectors, axisymmetric polynomial potentials',
          'formal_assistant_checked':False}
    cert['template_id']=hashlib.sha256(json.dumps(cert,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    return cert


def synthesize(direction,tolerance=F(1,10**6),max_leaves=64):
    dirs,w,_,_,g=assemble([direction],[[1]])
    return build(dirs,w,maximum(g,tolerance,max_leaves),'poisson_scalar_v29')


def verify(cert):
    try:return cert==build(cert['directions'],cert['weight'],cert['range_proof'],cert['method'])
    except (ValueError,KeyError,TypeError,IndexError,ZeroDivisionError):return False
