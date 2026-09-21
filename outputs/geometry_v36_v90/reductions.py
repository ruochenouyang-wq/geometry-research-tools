"""V24--V27: checked, value-preserving reductions with explicit witness lifts."""
import copy
from fractions import Fraction as F
import algebra as a
import compiler_v23 as c
import error_bounds as poly


def data(p):
    p=c.normalized(p);q=p['penalty']
    return list(map(F,p['q0'])),[list(map(F,v)) for v in p['directions']],F(q['constant']),list(map(F,q['linear'])),[list(map(F,r)) for r in q['hessian']]


def sphere_ground(p):
    p=c.normalized(p)
    if p['function_space']!='full_sphere':raise ValueError('Full sphere scope required')
    reduced=copy.deepcopy(p);reduced['function_space']='axisymmetric'
    return {'rule':'sphere_ground_v24','original':p,'reduced':reduced,'justification':'azimuthal_L2_envelope','witness_lift':'identity_axisymmetric_function'}


def absorb_and_eliminate(p):
    q0,dirs,cost,l,h=data(p);d=len(dirs);cost+=q0[0];q0[0]=F(0)
    for i in range(d):l[i]+=dirs[i][0];dirs[i][0]=F(0);dirs[i]=poly.trim(dirs[i])
    passive=[i for i,q in enumerate(dirs) if not any(q)];active=[i for i in range(d) if i not in passive]
    offset=[F(0)]*d;lift=[[F(i==j) for j in active] for i in range(d)]
    if passive:
        hh=[[h[i][j] for j in passive] for i in passive];ll=[l[i] for i in passive]
        if p['domain']['kind']=='all_real':
            if a.inertia(hh)!=[0,0,len(passive)]:raise ValueError('Eliminated block must be positive definite')
            inv=a.inverse(hh);base=[-v for v in a.matvec(inv,ll)]
            mixed=[[h[i][j] for j in active] for i in passive];slope=[[-v for v in row] for row in a.matmul(inv,mixed)]
            for z,i in enumerate(passive):offset[i]=base[z];lift[i]=slope[z]
            cost-=sum(x*y for x,y in zip(ll,a.matvec(inv,ll)))/2
            newl=[l[i]-sum(h[i][k]*v for k,v in zip(passive,a.matvec(inv,ll))) for i in active]
            newh=[[h[i][j]-sum(h[i][k]*inv[z][w]*h[m][j] for z,k in enumerate(passive) for w,m in enumerate(passive)) for j in active] for i in active]
            domain={'kind':'all_real'}
        else:
            if any(h[i][j] for i in passive for j in active):raise ValueError('Box elimination requires no active/passive coupling')
            box=[list(map(F,p['domain']['box'][i])) for i in passive];low,x=a.quadratic_min(0,ll,hh,box);cost+=low
            for i,v in zip(passive,x):offset[i]=v
            newl=[l[i] for i in active];newh=[[h[i][j] for j in active] for i in active]
            domain={'kind':'box','box':[p['domain']['box'][i] for i in active]}
    else:newl=l;newh=h;domain=p['domain']
    reduced=c.ir(q0,[dirs[i] for i in active],cost,newl,newh,domain,p['function_space'])
    if reduced==p:raise ValueError('No absorbable shift or passive coordinate')
    return {'rule':'passive_elimination_v25','original':p,'reduced':reduced,'active':active,'eliminated':passive,
            'lift_offset':list(map(str,offset)),'lift_matrix':[list(map(str,row)) for row in lift]}


def reflect(p):
    q0,dirs,cost,l,h=data(p);d=len(dirs)
    if p['domain']['kind']!='box' or any(q0[1::2]):raise ValueError('An even base potential and bounded box are required')
    signs=[]
    for q in dirs:
        even=any(q[::2]);odd=any(q[1::2])
        if even and odd:raise ValueError('Mixed parity direction')
        signs.append(-1 if odd else 1)
    odd=[i for i,s in enumerate(signs) if s<0]
    if not odd:raise ValueError('No reflection parameter')
    box=[list(map(F,r)) for r in p['domain']['box']]
    if any(box[i][0]!=-box[i][1] or l[i] for i in odd):raise ValueError('Reflection does not preserve domain or linear cost')
    if any(h[i][j] and signs[i]*signs[j]<0 for i in range(d) for j in range(d)):raise ValueError('Reflection does not preserve quadratic cost')
    axis=odd[0];box[axis][0]=F(0)
    reduced=c.ir(q0,dirs,cost,l,h,{'kind':'box','box':box},p['function_space'])
    return {'rule':'parameter_reflection_v26','original':p,'reduced':reduced,'signs':signs,'half_axis':axis,'witness_lift':'identity'}


def compress(p):
    q0,dirs,cost,l,h=data(p);d=len(dirs)
    if not d or p['domain']['kind']!='all_real':raise ValueError('Rank compression requires real parameters')
    if any(q[0] for q in dirs):raise ValueError('Absorb constant directions first')
    if a.inertia(h)!=[0,0,d]:raise ValueError('Positive definite penalty required')
    degree=max(map(len,dirs));matrix=[[q[k] if k<len(q) else F(0) for q in dirs] for k in range(1,degree)]
    indices=a.independent_rows(matrix);r=len(indices)
    if not 0<r<d:raise ValueError('No nontrivial rank reduction')
    C=[matrix[i] for i in indices];basisrows=[a.rectangular(a.transpose(C),row,r)[0] for row in matrix]
    newdirs=[poly.trim([F(0)]+[row[j] for row in basisrows]) for j in range(r)]
    inv=a.inverse(h);center=[-z for z in a.matvec(inv,l)];mean=a.matvec(C,center)
    effective=a.inverse(a.matmul(a.matmul(C,inv),a.transpose(C)))
    newl=[-z for z in a.matvec(effective,mean)];newcost=cost-sum(x*y for x,y in zip(l,a.matvec(inv,l)))/2+sum(x*y for x,y in zip(mean,a.matvec(effective,mean)))/2
    lift=a.matmul(a.matmul(inv,a.transpose(C)),effective);offset=[x-y for x,y in zip(center,a.matvec(lift,mean))]
    reduced=c.ir(q0,newdirs,newcost,newl,effective,{'kind':'all_real'},p['function_space'])
    return {'rule':'potential_rank_v27','original':p,'reduced':reduced,'coefficient_map':[list(map(str,row)) for row in C],
            'lift_offset':list(map(str,offset)),'lift_matrix':[list(map(str,row)) for row in lift]}


RULES={'sphere_ground_v24':sphere_ground,'passive_elimination_v25':absorb_and_eliminate,'parameter_reflection_v26':reflect,'potential_rank_v27':compress}


def verify(cert):
    try:return cert==RULES[cert['rule']](cert['original'])
    except (ValueError,KeyError,TypeError,IndexError,ZeroDivisionError):return False


def lift(cert,parameters):
    if not verify(cert):raise ValueError('Invalid reduction')
    if 'lift_matrix' not in cert:return list(map(F,parameters))
    return [F(x)+v for x,v in zip(cert['lift_offset'],a.matvec([list(map(F,row)) for row in cert['lift_matrix']],list(map(F,parameters))))]
