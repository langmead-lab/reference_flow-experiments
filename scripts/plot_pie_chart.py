import matplotlib.pyplot as plt

total = 0
num_chrs = 0
with open('num_highq_wgs.txt', 'r') as f:
    for line in f:
        if line.startswith('chr'):
            num_chrs += 1
        else:
            total += int(line.rstrip())

# Pie chart
labels = ['MAPQ >= 10', 'MAPQ < 10 and unaligned']
sizes = [total, num_chrs*200000-total]
#colors
# colors = ['orange', 'green']
# colors = ['#ff9999','#66b3ff','#99ff99','#ffcc99']

plt.rcParams.update({'font.size': 22})

fig1, ax1 = plt.subplots()
patches, texts, autotexts = ax1.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90)
# patches, texts, autotexts = ax1.pie(sizes, colors = colors, labels=labels, autopct='%1.1f%%', startangle=90)
# for text in texts:
#     text.set_color('grey')
# for autotext in autotexts:
#     autotext.set_color('grey')
# Equal aspect ratio ensures that pie is drawn as a circle
ax1.axis('equal')  
plt.tight_layout()
plt.show()

