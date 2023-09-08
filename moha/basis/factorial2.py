def fact2(n):
    if n ==1 or n==0 or n==-1:
        return 1
    else:
        return n*fact2(n-2)