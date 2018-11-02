cor = {42: 745621, 1: 27285, 40: 55352, 7: 2071, 30: 3386, 39: 12905, 22: 4525, 18: 1712, 6: 19175, 36: 3340, 17: 2066, 38: 6317, 24: 7403, 34: 2378, 27: 6012, 35: 3469, 31: 2753, 37: 6086, 25: 3308, 26: 7963, 12: 2673, 11: 1697, 33: 3022, 21: 2194, 16: 1159, 23: 2255, 32: 2612, 3: 114, 15: 1494, 14: 670, 0: 338, 4: 77, 2: 136, 5: 124, 8: 164}

incor = {42: 4, 1: 55245, 40: 0, 7: 105, 30: 8, 39: 0, 22: 8, 18: 4, 6: 819, 36: 0, 17: 18, 38: 0, 24: 0, 34: 0, 27: 1, 35: 0, 31: 3, 37: 0, 25: 4, 26: 6, 12: 35, 11: 35, 33: 0, 21: 6, 16: 2, 23: 2, 32: 1, 3: 31, 15: 17, 14: 4, 0: 1363, 4: 12, 2: 404, 5: 6, 8: 1}

if len(cor) != len(incor):
    print ("Error: diff length (%d/%d)" % (len(cor), len(incor)))

max_key = max(max(cor), max(incor))
min_key = min(min(cor), min(incor))

#for i in range(min_key, max_key + 1): # ascending
for i in range(max_key, min_key - 1, -1): # descending
    if cor.get(i):
        c = cor[i]
    else:
        c = 0
    if incor.get(i):
        ic = incor[i]
    else:
        ic = 0
    print (i, c, ic)

'''
scor = sorted(cor)
sincor = sorted(incor)
for i in range(len(cor)):
    print (scor[i], cor[scor[i]], incor[sincor[i]])
'''
