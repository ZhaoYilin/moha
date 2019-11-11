import numpy as np

def is_interaction(shape,listi,listj,operator,pbc):
    """

    """
    list_cell = list(operator.cell2-operator.cell1)
    bool2 = listi[3] == operator.site1
    bool3 = listj[3] == operator.site2
    if pbc==False:
        list_ij = list(np.array(listj[0:3])-np.array(listi[0:3]))
        bool1 = list_ij == list_cell
        return bool1 and bool2 and bool3
    elif pbc==True:
        list_ij_1 = list(np.array(listj[0:3])-np.array(listi[0:3])+np.array([shape[0],0,0]))
        list_ij_2 = list(np.array(listj[0:3])-np.array(listi[0:3])+np.array([0,shape[0],0]))
        list_ij_3 = list(np.array(listj[0:3])-np.array(listi[0:3])+np.array([0,0,shape[0]]))
        bool1 = list_ij_1==list_cell or list_ij_2==list_cell or list_ij_3==list_cell
        return bool1 and bool2 and bool3

def tensor_enumerate(tensor):
    n = -1
    result = []
    for i in range(tensor.shape[0]):
        for j in range(tensor.shape[1]):
            for k in range(tensor.shape[2]):
                for l in range(tensor.shape[3]):
                    n+=1
                    result.append((n,[i,j,k,l]))
    return result

def expectation(operator,site):
    state = site.states['occupied']
    if operator.__class__.__name__ == 'OneBodyTerm':
        if len(operator.terms) == 1:
            tmp = np.dot(site.operators[operator.terms[0]],state)
            tmp = np.dot(state,tmp)
        if len(operator.terms) == 2:
            tmp = np.dot(site.operators[operator.terms[1]],state)
            tmp = np.dot(site.operators[operator.terms[0]],tmp)
            tmp = np.dot(state,tmp)
    elif operator.__class__.__name__ == 'TwoBodyTerm':
        if len(operator.terms) == 2:
            tmp = np.dot(site.operators[operator.terms[1]],state)
            tmp = np.dot(site.operators[operator.terms[0]],tmp)
            tmp = np.dot(state,tmp)
        elif len(operator.terms) == 4:
            tmp1 = np.dot(site.operators[operator.terms[3]],state)
            tmp1 = np.dot(site.operators[operator.terms[2]],tmp1)
            tmp1 = np.dot(state,tmp1)
            tmp2 = np.dot(site.operators[operator.terms[1]],state)
            tmp2 = np.dot(site.operators[operator.terms[0]],tmp2)
            tmp2 = np.dot(state,tmp2)
            tmp = tmp1*tmp2
    return tmp*operator.parameter
