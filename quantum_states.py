import angmom as angmom
from itertools import groupby
from operator import itemgetter

def setup_NN_basis(j2min, j2max, tzmin, tzmax):
    basis = []
    for tz in range(tzmin,tzmax+1,1):
        for J in range(j2min,j2max+1,1):
            for S in range(0,2):
                for L in range(abs(J-S),J+S+1,1):
                    for T in range(abs(tz),2,1):
                        """
                        Looparna itererar med range(start,stop,steg)
                        från start till stop med ett "steg" i taget.
                        +1 pga. stop är inte inkluderat

                        """
                        if ((L+S+T)%2 != 0):
                            basis_state = {}
                            basis_state['tz'] = tz
                            basis_state['l']  = L
                            basis_state['pi'] = (-1)**L
                            basis_state['s']  = S
                            basis_state['j']  = J
                            basis_state['t']  = T
                            """
                            Om vilkoret (L+S+T)%2 != 0) uppfylls skapas en
                            "associative array/dictionary" basis_state
                            som använder kvanttalens "strings" som index och 
                            deras nuvarande nummer i loopen som värde.

                            """
                            print(basis_state)
                            basis.append(basis_state)
                            """
                            .append lägger till den skapade basen i slutet av basis listan

                            """


    print('len(basis) = ',len(basis))

    """
    Det som returneras är nu en lista med states.
     Basis har numrerad index, state är indexerad med "kvanttal".
    """
    return basis
 

def setup_NN_channels(basis):

    states = []
    
    for bra in basis:
        for ket in basis:
            """
            Denna syntax betyder bara "För varje element/item i listan basis"
            där elementen ges namnen "bra" och "ket", istället för numrerade index
            """
            
            if angmom.kroenecker_delta(bra,ket,'j','tz','s','pi'):
                """
                Här undersöks huruvida de relevanta storheterna är konserverade efter 
                spridning. Argumenten 'j','tz','s','pi' är återigen index för värdena
                på kvanttalen i tillstånds listorna "bra" och "ket", som representerar
                tillstånden innan och efter spridning. Kroeneckers delta returnerar
                true om kvanttalen är samma i "bra" och "ket".
                """
                
                state = {}
                
                state['l']  = bra['l']
                state['ll'] = ket['l']
                
                state['s']  = bra['s']
                state['j']  = bra['j']
                state['t']  = bra['t']
                state['tz'] = bra['tz']
                state['pi'] = bra['pi']
                states.append(state)

                """
                Nya "giltiga" tillstånd bildas av bastillstånden 
                och läggs till i dictionaryt states där 'kvanttal' är index
                och värdet är det nummer med motsvarande index.
                """


    """
    Short version: The following code sorts the states and put them into different 
    groups. The group it defined as the collection of states that have the same quantum numbers
    "s", "j", "tz", "pi". Each group is a channel and each state in the group is a block.

    https://stackoverflow.com/questions/773/how-do-i-use-itertools-groupby
    https://stackoverflow.com/questions/18595686/how-do-operator-itemgetter-and-sort-work

    """
    grouper = itemgetter("s", "j", "tz", "pi")
    NN_channels = []
    
    for key, grp in groupby(sorted(states, key = grouper), grouper):
        NN_channels.append(list(grp))
        

    for chn_idx, chn in enumerate(NN_channels):
        """
        enumerate lets you access both the key and value
        https://realpython.com/python-enumerate/

        """
        for block in chn:
            block.update({"chn_idx":chn_idx})
            """
            dictionary.update(iterable)
            Inserts the specified items to the dictionary,
            Where iterable objects has key-value pair.
            """
            
    for idx, chn in enumerate(NN_channels):
        print('channel %d: J=%d S=%d Tz=%d Pi=%d'%(idx,chn[0]['j'],chn[0]['s'],chn[0]['tz'],chn[0]['pi']))
        for block in chn:
            print('  block:',block)
            """
            What's the 0 in chn[0]["qn"]?

            Link to explaination of %d
            https://www.learnpython.org/en/String_Formatting
            """
    
    print('len(NN_channels) = ',len(NN_channels))
    return NN_channels
            
def lookup_NNchannel_idx(NN_channels, **kwargs):
    matching_indices = []
    for idx, chn in enumerate(NN_channels):
        for block in chn:
            if (kwargs.items() <= block.items()):
                matching_indices.append(idx)

    matching_indices = list(dict.fromkeys(matching_indices))

    return matching_indices
                
    
