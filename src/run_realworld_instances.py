import os

adolescent = "adolescent_health.lcip"
wikivote = "wiki-vote.lcip"
innovation = "innovation.lcip"
residense = "residense.lcip"
highschool = "highschool.lcip"

instances = list()

instances.append(highschool)
instances.append(residense)
instances.append(innovation)
instances.append(wikivote)
instances.append(adolescent)

print("network \talg \ttime \tnodes \tdualbound \tprimalbound \tgap")
for i in instances:
    cmd = '../bin/glcip -i ../in/realworld/' + i + ' -a arc -alpha .5 | grep -v "Academic license"'
    output = os.popen(cmd).read()
    print(i, "BC", output)
    cmd = '../bin/glcip -i ../in/realworld/' + i + ' -a arcwb -alpha .5 | grep -v "Academic license"'
    output = os.popen(cmd).read()
    print(i, "BC+", output)

print("alpha = .1")
for i in instances:
    cmd = '../bin/glcip -i ../in/realworld/' + i + ' -a arc -alpha .1 | grep -v "Academic license"'
    output = os.popen(cmd).read()
    print(i, "BC", output)
    cmd = '../bin/glcip -i ../in/realworld/' + i + ' -a arcwb -alpha .1 | grep -v "Academic license"'
    output = os.popen(cmd).read()
    print(i, "BC+", output)
