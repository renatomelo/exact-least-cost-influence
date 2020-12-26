import os

instances = list()
#inserted in a non-decreasing order by the number of vertices
instances.append("highschool.lcip")
instances.append("residense.lcip")
instances.append("innovation.lcip")
instances.append("wiki-vote.lcip")
instances.append("adolescent_health.lcip")
instances.append("advogato.lcip")
instances.append("dblp.lcip")
instances.append("cora.lcip")
instances.append("HepTh.lcip")

#paremeters to vary
alphas = list()

alphas.append("1")
alphas.append(".5")
alphas.append(".1")

for alpha in alphas: 
    print("alpha = ", alpha)
    print("network \talg \ttime \tnodes \tdualbound \tprimalbound \tgap")
    for i in instances:
        cmd = '../bin/glcip -i ../data/realworld/' + i + ' -a arc -alpha '+ alpha +' | grep -v "Academic license"'
        output = os.popen(cmd).read()
        output = output.replace("\n", "")
        print(i, "BC \t", output)
        cmd = '../bin/glcip -i ../data/realworld/' + i + ' -a arcwb -alpha '+ alpha +' | grep -v "Academic license"'
        output = os.popen(cmd).read()
        output = output.replace("\n", "")
        print(i, "BC+\t", output)

""" print("alpha = .1")
for i in instances:
    cmd = '../bin/glcip -i ../data/realworld/' + i + ' -a arc -alpha .1 | grep -v "Academic license"'
    output = os.popen(cmd).read()
    print(i, "BC", output)
    cmd = '../bin/glcip -i ../data/realworld/' + i + ' -a arcwb -alpha .1 | grep -v "Academic license"'
    output = os.popen(cmd).read()
    print(i, "BC+", output) """
