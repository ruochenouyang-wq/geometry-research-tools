"""V23: exact polynomial-expression compiler with explicit geometric scope.

This accepts a documented expression language, not arbitrary natural language.
"""
import ast
import copy
import re
from fractions import Fraction as F
import spectral_certifier as base


def add(a,b):
    out=dict(a)
    for e,v in b.items():
        out[e]=out.get(e,F(0))+v
        if not out[e]:del out[e]
    return out


def mul(a,b):
    out={}
    for e,v in a.items():
        for f,w in b.items():
            k=tuple(x+y for x,y in zip(e,f))
            if sum(k)>16:raise ValueError('Expression degree budget')
            out[k]=out.get(k,F(0))+v*w
    if len(out)>2048:raise ValueError('Expression term budget')
    return {e:v for e,v in out.items() if v}


def polynomial(text,names):
    if not isinstance(text,str) or len(text)>8192:raise ValueError('Expression size')
    zero=(0,)*len(names);tree=ast.parse(text,mode='eval')
    if sum(1 for _ in ast.walk(tree))>512:raise ValueError('Expression tree budget')
    def visit(node):
        if isinstance(node,ast.Constant) and type(node.value) is int:
            if abs(node.value)>10**12:raise ValueError('Integer input budget')
            return {zero:F(node.value)} if node.value else {}
        if isinstance(node,ast.Name) and node.id in names:
            e=list(zero);e[names.index(node.id)]=1;return {tuple(e):F(1)}
        if isinstance(node,ast.UnaryOp) and isinstance(node.op,(ast.UAdd,ast.USub)):
            return {e:v*(-1 if isinstance(node.op,ast.USub) else 1) for e,v in visit(node.operand).items()}
        if isinstance(node,ast.BinOp):
            a=visit(node.left)
            if isinstance(node.op,ast.Pow):
                if not isinstance(node.right,ast.Constant) or type(node.right.value) is not int or not 0<=node.right.value<=12:raise ValueError('Power must be an integer 0..12')
                out={zero:F(1)}
                for _ in range(node.right.value):out=mul(out,a)
                return out
            b=visit(node.right)
            if isinstance(node.op,ast.Add):return add(a,b)
            if isinstance(node.op,ast.Sub):return add(a,{e:-v for e,v in b.items()})
            if isinstance(node.op,ast.Mult):return mul(a,b)
            if isinstance(node.op,ast.Div) and set(b)=={zero} and b[zero]:return {e:v/b[zero] for e,v in a.items()}
        raise ValueError('Unsupported expression; use integers, rational constants, declared names, + - * / and **')
    return visit(tree.body)


def ir(q0,directions,c,l,h,domain,scope='axisymmetric'):
    d=len(directions)
    if not 0<=d<=4 or len(l)!=d or len(h)!=d or any(len(row)!=d for row in h):raise ValueError('Parameter dimensions')
    q0=base.potential(list(map(str,q0)));directions=[base.potential(list(map(str,p))) for p in directions]
    h=[list(map(F,row)) for row in h];l=list(map(F,l))
    if any(h[i][j]!=h[j][i] for i in range(d) for j in range(d)):raise ValueError('Non-symmetric Hessian')
    if scope not in ('axisymmetric','full_sphere'):raise ValueError('Unsupported function space')
    if domain=={'kind':'all_real'}:domain=dict(domain)
    elif set(domain)=={'kind','box'} and domain['kind']=='box':
        box=[list(map(F,row)) for row in domain['box']]
        if len(box)!=d or any(len(row)!=2 or row[0]>=row[1] for row in box):raise ValueError('Nonempty parameter box required')
        domain={'kind':'box','box':[list(map(str,row)) for row in box]}
    else:raise ValueError('Unsupported parameter domain')
    return {'geometry':'unit_sphere','function_space':scope,'eigenvalue_index':1,'q0':list(map(str,q0)),
            'directions':[list(map(str,p)) for p in directions],
            'penalty':{'constant':str(F(c)),'linear':list(map(str,l)),'hessian':[list(map(str,row)) for row in h]},'domain':domain}


def normalized(p):
    if set(p)!={'geometry','function_space','eigenvalue_index','q0','directions','penalty','domain'} or p['geometry']!='unit_sphere' or type(p['eigenvalue_index']) is not int or p['eigenvalue_index']!=1:raise ValueError('Only ground energy of the unit sphere is supported')
    q=p['penalty'];expected=ir(p['q0'],p['directions'],q['constant'],q['linear'],q['hessian'],p['domain'],p['function_space'])
    if expected!=p:raise ValueError('IR must be canonical')
    return expected


def compile_goal(raw):
    expected={'geometry','function_space','eigenvalue_index','parameters','potential','penalty','domain','threshold'}
    if not isinstance(raw,dict) or set(raw)!=expected:raise ValueError('Explicit goal fields required')
    if raw['geometry']!='unit_sphere' or type(raw['eigenvalue_index']) is not int or raw['eigenvalue_index']!=1:raise ValueError('Unsupported geometry or spectral index')
    names=raw['parameters'];d=len(names)
    if not isinstance(names,list) or not 0<=d<=4 or len(set(names))!=d or any(not isinstance(n,str) or re.fullmatch('[a-zA-Z][a-zA-Z0-9_]*',n) is None or n=='t' for n in names):raise ValueError('Distinct parameter names required')
    potential=polynomial(raw['potential'],names+['t']);q0=[F(0)]*7;directions=[[F(0)]*7 for _ in names]
    for e,v in potential.items():
        if sum(e[:-1])>1 or e[-1]>6:raise ValueError('Potential must be affine in parameters and degree <= 6 in t')
        if not sum(e[:-1]):q0[e[-1]]+=v
        else:directions[e[:-1].index(1)][e[-1]]+=v
    cost=polynomial(raw['penalty'],names);c=F(0);l=[F(0)]*d;h=[[F(0)]*d for _ in names]
    for e,v in cost.items():
        if sum(e)>2:raise ValueError('Penalty must be quadratic')
        if sum(e)==0:c+=v
        elif sum(e)==1:l[e.index(1)]+=v
        elif 2 in e:h[e.index(2)][e.index(2)]+=2*v
        else:
            i,j=[i for i,z in enumerate(e) if z];h[i][j]+=v;h[j][i]+=v
    problem=ir(q0,directions,c,l,h,raw['domain'],raw['function_space'])
    bound=F(str(raw['threshold']))
    return {'method':'expression_compiler_v23','source':copy.deepcopy(raw),'problem':problem,'threshold':str(bound),'relation':'>='}


def verify(cert):
    try:return cert==compile_goal(cert['source'])
    except (ValueError,KeyError,TypeError,SyntaxError,IndexError,ZeroDivisionError,RecursionError):return False
