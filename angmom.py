# Checks triangular relation
# |a-b| <= ab <= a+b TRUE
# ELSE: FALSE
def triag(a, b, ab):
    if( ab < abs(a - b) ):
        return False
    if ( ab > a + b ):
        return False
    return True

def kroenecker_delta(bra,ket,*args):
    """
    *args används när antalet argument är
    okänt, kan varieras vid runtime
    """
    for ar in args:
        if bra[ar] != ket[ar]:
            return False
    return True
