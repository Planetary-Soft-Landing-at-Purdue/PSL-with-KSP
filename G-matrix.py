G_txt = ''
for r in range(len(G)):
    for c in G[r]:
        if c == 0:
            G_txt += ' '*7 + '-'
        else:
            G_txt += ' '*(8-len(str(round(-c, 3)))) + str(round(-c, 3))
    G_txt += '  |' + ' '*(8-len(str(round(h[r], 3)))) + str(round(h[r], 3)) +  '\n'
txtFile = open("G_matrix.txt", 'w')
txtFile.write(G_txt)
txtFile.close()