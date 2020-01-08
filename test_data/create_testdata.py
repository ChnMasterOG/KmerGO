import random
int2word = ['A', 'C', 'G', 'T']
reads_length = 80
reads_number = 1000
specific_seq_number = 10
specific_seq_length = 150
specific_seq_list = ['' for i in range(specific_seq_number)]

# create random specific_seq
for i in range(specific_seq_number):
    for j in range(specific_seq_length):
        specific_seq_list[i] += int2word[random.randint(0, 3)]

print(specific_seq_list)

for z in range(25): # Group A
    # create reads
    reads = []
    for i in range(specific_seq_number):
        for j in range(int(specific_seq_length / reads_length) + 2):
            reads.append(specific_seq_list[i][min(j*reads_length, specific_seq_length-reads_length):min((j+1)*reads_length, specific_seq_length)])

    while len(reads) < reads_number:
        s = ''
        for i in range(reads_length):
            s += int2word[random.randint(0, 3)]
        reads.append(s)

    random.shuffle(reads)

    f = open('group_A/A' + str(z+1) + '.fasta', 'w')
    for i in range(len(reads)):
        f.write('>' + str(i+1) + '\n')
        f.write(reads[i] + '\n')
    f.close()

for z in range(25): # Group B
    # create reads
    reads = []
    
    while len(reads) < reads_number:
        s = ''
        for i in range(reads_length):
            s += int2word[random.randint(0, 3)]
        reads.append(s)

    f = open('group_B/B' + str(z+1) + '.fasta', 'w')
    for i in range(len(reads)):
        f.write('>' + str(i+1) + '\n')
        f.write(reads[i] + '\n')
    f.close()
