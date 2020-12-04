import os

instances = list()

#read the list_of_instance.txt
input_file = open("../data/list_of_instances.txt", "r")

#iterate over the lines of a file
for line in input_file:
    #add to the list the file's name
    instances.append(line.replace("\n",''))

input_file.close()

#paremeters to vary
alphas = list()
algorithms = list()

alphas.append(1)
alphas.append(.5)
alphas.append(.1)

algorithms.append("arc")
algorithms.append("arcwb")
#algorithms.append("covcg")
#algorithms.append("cov")

for alpha in alphas: 
    print("alpha = ", alpha)
    i = 0
    print("network \talg \ttime \tnodes \tdualbound \tprimalbound \tgap")
    while i < len(instances):
        for a in algorithms:
            values = list([0] * 6)
            for j in range(5):
                #print(a, instances[i + j])
                name = instances[i + j]
                cmd = '../bin/glcip -i ../data/socnet-instances/' + name + ' -a '+ a + ' -alpha ' + str(alpha) + ' | grep -v "Academic license"'
                #print(cmd)
                output = os.popen(cmd).read()
                output = output.replace("\n", "")
                splited = output.split("\t")
                for x in range(len(values)):
                    values[x] = values[x] + float(splited[x])
                print(name, a, output)
            #compute the average here
            print("Average\t\t", end = ' ')
            for x in values:
                print("\t", x/5, end = '')
            print()
        i = i + j + 1
