"""V31: transport a checked family lemma to affine amplitudes and dominated potentials."""
from fractions import Fraction as F
import compiler_v23 as compiler
import reductions
import algebra as a
import family_v29 as family
import error_bounds as e
import v05_range


def apply(p,template,remainder_proof=None):
    p=compiler.normalized(p)
    if not family.verify(template):raise ValueError('Invalid reusable lemma')
    q0,dirs,c,l,h=reductions.data(p);basis=[list(map(F,q)) for q in template['directions']];r=len(basis);d=len(dirs)
    degree=max([len(q0)]+[len(q) for q in basis+dirs])
    matrix=[[q[k] if k<len(q) else F(0) for q in basis] for k in range(1,degree)]
    if len(a.independent_rows(matrix))!=r:raise ValueError('Template directions must be independent modulo constants')
    pivots=a.independent_rows(matrix)
    def coordinates(q,exact):
        target=[q[k] if k<len(q) else F(0) for k in range(1,degree)]
        if exact:return a.rectangular(matrix,target,r)[0]
        return a.solve([matrix[k] for k in pivots],[target[k] for k in pivots])
    offset=coordinates(q0,False);columns=[coordinates(q,True) for q in dirs];mapping=[[col[j] for col in columns] for j in range(r)]
    remainder=q0[:]
    for z,q in zip(offset,basis):remainder=e.add(remainder,e.scale(q,-z))
    if remainder_proof is None:remainder_proof=v05_range.make_range(remainder)
    remlo,_=v05_range.verify_bounds(remainder,remainder_proof)
    constant_dirs=[q[0]-sum(columns[i][j]*basis[j][0] for j in range(r)) for i,q in enumerate(dirs)]
    means=list(map(F,template['means']));g=[list(map(F,row)) for row in template['quadratic_bound']];gz=a.matvec(g,offset)
    newc=c+remlo+sum(v*w for v,w in zip(means,offset))-sum(v*w for v,w in zip(offset,gz))
    newl=[l[i]+constant_dirs[i]+sum(mapping[j][i]*(means[j]-2*gz[j]) for j in range(r)) for i in range(d)]
    newh=[[h[i][j]-2*sum(mapping[z][i]*g[z][w]*mapping[w][j] for z in range(r) for w in range(r)) for j in range(d)] for i in range(d)]
    box=[list(map(F,row)) for row in p['domain']['box']] if p['domain']['kind']=='box' else None
    low,where=a.quadratic_min(newc,newl,newh,box)
    return {'method':'lemma_transport_v31','problem':p,'template':template,'amplitude_offset':list(map(str,offset)),
            'amplitude_matrix':[list(map(str,row)) for row in mapping],'remainder_polynomial':list(map(str,remainder)),
            'remainder_proof':remainder_proof,'derived_quadratic':{'constant':str(newc),'linear':list(map(str,newl)),'hessian':[list(map(str,row)) for row in newh]},
            'lower_bound':str(low),'relaxation_minimizer':list(map(str,where)),
            'formal_assistant_checked':False}


def verify(cert):
    try:return cert==apply(cert['problem'],cert['template'],cert['remainder_proof'])
    except (ValueError,KeyError,TypeError,IndexError,ZeroDivisionError):return False
